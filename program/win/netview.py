#!/usr/bin/env python
# NetView P v.0.7 - Windows
# Dependencies: PLINK
# Eike Steinig
# Zenger Lab, JCU
# https://github.com/esteinig/netview

# Multiprocessing module for Windows

import os
import time
import shutil
import argparse
import subprocess
import numpy as np
import multiprocessing as mp
import scipy.sparse.csgraph as csg
import scipy.spatial.distance as sd

from sklearn.neighbors import NearestNeighbors

def main():

    commands = CommandLine()

    dat = Data()
    dat.prefix = commands.arg_dict['prefix']
    dat.ploidy = commands.arg_dict['ploidy']
    dat.missing = commands.arg_dict['missing']

    if commands.arg_dict['visual']:
        print('\nGenerated node attribute files only.\n')
        dat.readData(commands.arg_dict['attribute_file'], f='attributes', sep=',')
        dat.writeData(f='attributes')
        makeProject(commands.arg_dict['project'] + '_attributes', commands.arg_dict['prefix'])
        exit(1)

    print()
    print(get_time() + "\t" + "---------------------------------")
    print(get_time() + "\t" + "         NETVIEW P v.0.7         ")
    print(get_time() + "\t" + "---------------------------------")
    print(get_time() + "\t" + "File =", commands.arg_dict['data_file'].upper())

    if commands.arg_dict['plink']:
        dat.filetype = 'plink'
        dat.readData(commands.arg_dict['data_file'], f='plink', sep=commands.arg_dict['sep'])
    elif commands.arg_dict['snps']:
        dat.filetype = 'snps'
        dat.readData(commands.arg_dict['data_file'], f='snp_mat', sep=commands.arg_dict['sep'])
    else:
        dat.filetype = 'dist'
        dat.readData(commands.arg_dict['data_file'], f='matrix', sep=commands.arg_dict['sep'])

    dat.readData(commands.arg_dict['attribute_file'], f='attributes', sep=',')

    if dat.ploidy == 'diploid':
        nsnp = dat.nSNP//2
    else:
        nsnp = dat.nSNP

    print(get_time() + "\t" + "N =", str(dat.n).upper())
    print(get_time() + "\t" + "SNPs =", str(nsnp).upper())
    print(get_time() + "\t" + "Ploidy =", dat.ploidy.upper())
    print(get_time() + "\t" + "---------------------------------")
    print(get_time() + "\t" + "Quality Control =", str(commands.arg_dict['qc']).upper())

    pipeline = Analysis(dat)

    qc = False
    if commands.arg_dict['qc'] and pipeline.data.filetype != 'dist':

        qc_params = {'--mind': commands.arg_dict['mind'],
                     '--geno': commands.arg_dict['geno'],
                     '--maf': commands.arg_dict['maf'],
                     '--hwe': commands.arg_dict['hwe']}

        pipeline.runPLINK(qc_parameters=qc_params, quality=True)
        qc = True

    if commands.arg_dict['mat'] and pipeline.data.filetype != 'dist':
        pipeline.getDistance(distance=commands.arg_dict['distance'])
        pipeline.data.writeData(file=commands.arg_dict['prefix'] + '_mat.dist', f='matrix')
        makeProject(commands.arg_dict['project'] + '_dist', commands.arg_dict['prefix'])
        print(get_time() + "\t" + "---------------------------------\n")
        exit(1)
    elif commands.arg_dict['mat'] and pipeline.data.filetype == 'dist':
        print('\nError. Input is already a Distance Matrix.\n')
        exit(1)

    if not commands.arg_dict['off']:
        if pipeline.data.filetype != 'dist':
            pipeline.getDistance(distance=commands.arg_dict['distance'])
        pipeline.runNetView(tree=commands.arg_dict['tree'], start=commands.arg_dict['start'],
                            stop=commands.arg_dict['stop'], step=commands.arg_dict['step'],
                            algorithm=commands.arg_dict['algorithm'])

    if qc:
        pipeline.updateNodeAttributes(commands.arg_dict['attribute_file'])
    pipeline.data.writeData(f='attributes')

    makeProject(commands.arg_dict['project'], commands.arg_dict['prefix'])

    print(get_time() + "\t" + "---------------------------------\n")

def makeProject(project, prefix):

    cwd = os.getcwd()

    project_path = os.path.realpath(os.path.join(os.getcwd(), project))
    plink_path = os.path.realpath(os.path.join(project_path, 'plink'))
    network_path = os.path.realpath(os.path.join(project_path, 'networks'))
    other_path = os.path.realpath(os.path.join(project_path, 'other'))
    node_path = os.path.realpath(os.path.join(project_path, 'nodes'))

    if os.path.exists(project_path):
        shutil.rmtree(project_path)

    architecture = [project_path, plink_path, network_path, other_path, node_path]

    for directory in architecture:
        try:
            os.makedirs(directory)
        except OSError:
            if not os.path.isdir(directory):
                raise

    for name in os.listdir(cwd):
        if name.endswith('.edges'):
            pathname = os.path.join(cwd, name)
            if os.path.isfile(pathname):
                shutil.move(pathname, network_path)
        if name.endswith('.dist'):
            pathname = os.path.join(cwd, name)
            if os.path.isfile(pathname):
                shutil.move(pathname, other_path)
        if name.endswith('.nat'):
            pathname = os.path.join(cwd, name)
            if os.path.isfile(pathname):
                shutil.move(pathname, node_path)
        elif name.startswith(prefix + '_plink_in'):
            pathname = os.path.join(cwd, name)
            if os.path.isfile(pathname):
                os.remove(pathname)
        elif name.startswith(prefix + '_plink'):
            pathname = os.path.join(cwd, name)
            if os.path.isfile(pathname):
                shutil.move(pathname, plink_path)
        elif name.endswith('_qc.csv'):
            pathname = os.path.join(cwd, name)
            if os.path.isfile(pathname):
                shutil.move(pathname, other_path)

#### Functions for Multiprocessing ####

def netview(matrix, k, mst, algorithm, tree):

    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm=algorithm).fit(matrix)
    adj_knn = nbrs.kneighbors_graph(matrix).toarray()
    np.fill_diagonal(adj_knn, 0)
    adj_mknn = (adj_knn == adj_knn.T) * adj_knn

    if tree:
        adj = mst + adj_mknn
    else:
        adj = adj_mknn

    adjacency = np.tril(adj)
    mst_edges = np.argwhere(adjacency < 1)
    adjacency[adjacency > 0] = 1.
    edges = np.argwhere(adjacency != 0)
    weights = matrix[edges[:, 0], edges[:, 1]]

    return [k, edges, weights, adjacency, mst_edges]

def netview_callback(k):
    print(get_time() + "\t" + '              k=' + str(k[0]))

def get_time():
    return time.strftime("[%H:%M:%S]")

#### Command Line Module ####

class CommandLine:

    def __init__(self):

        self.parser = argparse.ArgumentParser(description='NetView P v0.7', add_help=True)
        self.setParser()
        self.args = self.parser.parse_args()

        self.arg_dict = vars(self.args)

    def setParser(self):

        data_type = self.parser.add_mutually_exclusive_group(required=True)

        # Required Options

        self.parser.add_argument('-f', dest='data_file', required=True, type=str,
                                help="Name of Data File")
        data_type.add_argument('-p', dest='plink', action='store_true',
                               help="PLINK format (.ped/.map)")
        data_type.add_argument('-s', dest='snps', action='store_true',
                               help="SNP matrix (N x SNPs)")
        data_type.add_argument('-m', dest='dist', action='store_true',
                               help="Distance matrix (N x N)")


        self.parser.add_argument('-a', dest='attribute_file', default='', type=str, required=True,
                                 help="Node attribute file (.csv)")

        # MAIN Options

        self.parser.add_argument('--quality', dest='qc', action='store_true', default=False,
                                 help="Quality control in PLINK (OFF)")
        self.parser.add_argument('--distance', dest='distance', default='asd', type=str,
                               help="Distance measure for SNPs: hamming, asd, correlation... (asd)")
        self.parser.add_argument('--algorithm', dest='algorithm', default='auto', type=str,
                               help="Algorithm for NN: auto, ball_tree, kd_tree, brute (brute)")
        self.parser.add_argument('--mst-off', dest='tree', action='store_false', default=True,
                                 help="Disable minimum spanning tree (OFF)")

        self.parser.add_argument('--ploidy', dest='ploidy', default='diploid', type=str,
                                 help="Set ploidy: haploid, diploid (diploid.")
        self.parser.add_argument('--missing', dest='missing', default='0', type=str,
                                 help="Set missing character (0)")

        self.parser.add_argument('--prefix', dest='prefix', default='project', type=str,
                                 help="Set prefix (project)")
        self.parser.add_argument('--project', dest='project', default=time.strftime("%d-%m-%Y_%H-%M-%S"), type=str,
                                 help="Output project name (timestamp)")
        self.parser.add_argument('--sep', dest='sep', default='\t', type=str,
                                 help="Delimiter for data file (\\t).")

        # PARAMETER Options

        self.parser.add_argument('--mind', dest='mind', default=0.1, type=float,
                                 help="Filter samples > missing rate (0.1)")
        self.parser.add_argument('--geno', dest='geno', default=0.1, type=float,
                                 help="Filter SNPs > missing rate (0.1)")
        self.parser.add_argument('--maf', dest='maf', default=0.01, type=float,
                                 help="Filter SNPs < minor allele frequency (0.01)")
        self.parser.add_argument('--hwe', dest='hwe', default=0.001, type=float,
                                 help="Filter SNPs failing HWE test at P < (0.001)")
        self.parser.add_argument('--start', dest='start', default=10, type=int,
                                 help="Start at k = (10)")
        self.parser.add_argument('--stop', dest='stop', default=40, type=int,
                                 help="Stop at k = (40)")
        self.parser.add_argument('--step', dest='step', default=10, type=int,
                                 help="Step by k = (10)")

        # PIPELINE Options

        self.parser.add_argument('--visual', dest='visual', action='store_true', default=False,
                                 help="Node attributes ONLY (OFF)")
        self.parser.add_argument('--off', dest='off', action='store_true', default=False,
                                 help="Switch off NetView and run only QC (OFF).")
        self.parser.add_argument('--matrix', dest='mat', action='store_true', default=False,
                                 help="Generate distance matrix ONLY (OFF).")

#### Data Module ####

class Data:
    
    ### DATA ATTRIBUTES ###
    
    def __init__(self):
        
        self.prefix = "project"
        self.ploidy = 'diploid'
        self.missing = "0"
        
        self.n = 0
        self.nSNP = 0
        self.ids = []  # IDs

        self.alleles = []
        self.snps = np.arange(5)  # Array (N x SNPs)
        self.biodata = []  # List/Alignment of BioPython SeqRecords
        
        self.meta_data = {}
        self.snp_data = {}
        
        self.matrices = {}
        self.networks = {}

        self.matrix = np.arange(5)  # Current Matrix
        
        self.netview_runs = 0
        self.filetype = ''

    ### DATA READER ###
    
    def readData(self, file, f, sep="\t", header=False, add_col=0):

        def _read_nexus(file, sep=sep):
                       
            snp_position = []
            snps = []
            matrix = False
            for line in file:
                content = line.strip().split(sep)
                if matrix == True:
                    if ";" in line:
                        break
                    snp_position.append(content[0])
                    snps.append(content[1:])
                else:
                    if "dimensions" in line:
                        self.n = int(content[1].split("=")[1])
                        self.nSNP = int(content[2].split("=")[1][:-1])
                    elif "taxlabels" in line:
                        self.ids = content[1:]
                    elif "matrix" in line:
                        matrix = True

            self.snps = np.array([list(i) for i in zip(*snps)]) # ordered by N
            self.snp_data['snp_id'] = [''.join(p.split("_")[:-1]) for p in snp_position]
            self.snp_data['snp_position'] = [p.split("_")[-1] for p in snp_position]
            self.filetype = 'nexus'

        def _read_raxml(file, sep=sep):

            header = []
            ids = []
            snps = []
            for line in file:
                content = line.strip().split(sep)
                if header:
                    ids.append(content[0])
                    snps.append(content[1])
                else:
                    header = content
                    self.n = header[0]
                    self.nSNP = header[1]
            
            self.ids = ids
            self.snps = np.array(snps)
            self.filetype = 'raxml'

        def _read_plink(file, filename, sep=sep):

            map_name = filename.split(".")[0] + ".map"
            map_file = open(map_name)

            ids = []
            meta = []
            snps = []
            for line in file:
                content = line.strip().split(sep)
                ids.append(content[1])
                snps.append(content[6:])
                meta.append(content[:6])
            
            self.ids = ids
            self.snps = np.array(snps)
            self.nSNP = len(self.snps[0])
            self.n = len(self.ids)
        
            self.meta_data["pop"] = [i[0] for i in meta]
            self.meta_data["dam"] = [i[2] for i in meta]
            self.meta_data["sire"] = [i[3] for i in meta]
            self.meta_data["sex"] = [i[4] for i in meta]
            self.meta_data["phenotype"] = [i[5] for i in meta]
                        
            map_content = [line.strip().split() for line in map_file]
            map_content = list(zip(*map_content))
            self.snp_data['snp_chromosome'] = list(map_content[0])
            self.snp_data['snp_id'] = list(map_content[1])
            self.snp_data['snp_genetic_distance'] = list(map_content[2])
            self.snp_data['snp_position'] = list(map_content[3])

            map_file.close()

            self.filetype = 'plink'

        def _read_matrix(file, header=header, add_col=add_col, sep=sep):
            
            content = [line.strip().split(sep)[add_col:] for line in file]
            if header:
                content = content[1:]

            matrix = np.array([list(map(float, ind)) for ind in content])
            
            self.matrix = matrix
            self.matrices['input'] = matrix

            return matrix

        def _read_snp_mat(file, sep):

            matrix = np.array([line.strip().split(sep) for line in file])
            self.snps = matrix
            self.n = len(matrix[:, 1])
            self.nSNP = len(matrix[1, :])
            if self.ploidy == 'diploid':
                self.snp_data['snp_id'] = [str(i) for i in range(self.nSNP//2)]
            else:
                self.snp_data['snp_id'] = [str(i) for i in range(self.nSNP)]

        def _read_attributes(file, sep=sep):
            content = [line.strip().split(sep) for line in file]
            head = content[0]
            content = list(zip(*content[1:]))
            for i in range(len(head)):
                self.meta_data[head[i]] = content[i]
            self.ids = list(content[0])

        ## Main Read ##
        
        infile = open(file)
        f = f.lower()
        
        if f == "nexus":
            _read_nexus(infile, sep)
        elif f =="raxml":
            _read_raxml(infile, sep)
        elif f == "plink":
            _read_plink(infile, file, sep)
        elif f == "matrix":
            matrix = _read_matrix(infile, header, add_col, sep)
        elif f == 'snp_mat':
            _read_snp_mat(infile, sep)
        elif f == 'attributes':
            _read_attributes(infile, sep)
        else:
            print("File format not supported.")
            raise IOError
        
        infile.close()

        if f != 'attributes':
            alleles = np.unique(self.snps).tolist()
            if self.missing in alleles:
                alleles.remove(self.missing)
            self.alleles = alleles

        if f == 'matrix':
            return matrix

    ### DATA WRITER ###
        
    def writeData(self, f, file='data.out', sep="\t"):

        def _write_raxml(outfile, sep):
            
            outfile.write(str(self.n) + sep + str(self.nSNP) + "\n")
            for i in range(self.n):
                outfile.write(self.ids[i] + sep + ''.join(self.snps[i]) + "\n")
        
        def _write_nexus(outfile, sep):
                
            taxlabels = " ".join(self.ids)
            header = '#nexus\nbegin data;\ndimensions ntax=' + str(self.n) + ' nchar=' + str(self.nSNP) + \
            ';\nformat symbols="AGCT" gap=. datatype=nucleotide;\ntaxlabels ' + taxlabels + ';\nmatrix\n'
            tail = ";\nend;"

            snps = list(zip(*self.snps))

            outfile.write(header)
            for i in range(self.nSNP):
                if 'snp_chromosome' in self.snp_data.keys():
                    outfile.write(self.snp_data['snp_chromosome'][i] + "_")
                else:
                    outfile.write(sep)
                if 'snp_id' in self.snp_data.keys():
                    outfile.write(self.snp_data['snp_id'][i] + sep)
                else:
                    outfile.write("SNP" + str(i) + sep)
                
                outfile.write(sep.join(snps[i]) + "\n")

            outfile.write(tail)

        def _write_plink(outfile, filename, sep):
            
            mapname = filename.split('.')[0] + ".map"

            for i in range(self.n):
                if 'pop' in self.meta_data.keys():
                    outfile.write(self.meta_data['pop'][i] + sep)
                else:
                    outfile.write("NA" + sep)
                if self.ids:
                    outfile.write(self.ids[i] + sep)
                else:
                    outfile.write("N" + str(i+1) + sep)
                if 'dam' in self.meta_data.keys():
                    outfile.write(self.meta_data['dam'][i] + sep)
                else:
                    outfile.write("0" + sep)
                if 'sire' in self.meta_data.keys():
                    outfile.write(self.meta_data['sire'][i] + sep)
                else:
                    outfile.write("0" + sep)
                if 'sex' in self.meta_data.keys():
                    outfile.write(self.meta_data['sex'][i] + sep)
                else:
                    outfile.write("0" + sep)
                if 'phenotype' in self.meta_data.keys():
                    outfile.write(self.meta_data['phenotype'][i] + sep)
                else:
                    outfile.write("0" + sep)
                    
                outfile.write(sep.join(self.snps[i]) + "\n")
            
            map_file = open(mapname, "w")

            if 'snp_id' in self.snp_data:
                for i in range(len(self.snp_data['snp_id'])):
                    if 'snp_chromosome' in self.snp_data.keys():
                        map_file.write(self.snp_data['snp_chromosome'][i] + sep)
                    else:
                        map_file.write("0" + sep)
                    if 'snp_id' in self.snp_data.keys():
                        map_file.write(self.snp_data['snp_id'][i] + sep)
                    else:
                        map_file.write("SNP" + str(i+1) + sep)
                    if 'snp_genetic_distance' in self.snp_data.keys():
                        map_file.write(self.snp_data['snp_genetic_distance'][i] + sep)
                    else:
                        map_file.write("0" + sep)
                    if 'snp_position' in self.snp_data.keys():
                        map_file.write(self.snp_data['snp_position'][i] + sep + "\n")
                    else:
                        map_file.write("0" + sep + "\n")

            map_file.close()
        
        def _write_metadata(outfile, sep):
        
            outfile.write("#" + sep + "n=" + str(self.n) + sep + "nSNP=" +
                          str(self.nSNP) + sep + "(" + self.ploidy + ")\n")
            ordered_keys = sorted([key for key in self.meta_data.keys()])
            outfile.write("Isolate")
            for key in ordered_keys:
                outfile.write(sep + key)
            outfile.write("\n")
            for i in range(self.n):
                if self.ids:
                    outfile.write(self.ids[i])
                else:
                    outfile.write("N" + str(1))
                for key in ordered_keys:
                    outfile.write(sep + self.meta_data[key][i])
                outfile.write("\n")
        
        def _write_snpdata(outfile, sep):
            
            outfile.write("#" + sep + "n=" + str(self.n) + sep + "nSNP=" +
                          str(self.nSNP) + sep + "(" + self.ploidy + ")\n")
            
            snp_data = dict(self.snp_data)
            ordered_keys = sorted([key for key in snp_data.keys()])
            
            outfile.write("SNP" + sep)
            for key in ordered_keys:
                outfile.write(sep + key)
            outfile.write("\n")
            
            for i in range(self.nSNP):
                outfile.write("SNP_" + str(i))
                for key in ordered_keys:
                    outfile.write(sep + snp_data[key][i])
                outfile.write("\n")

        def _write_attributes():

            for key, value in self.meta_data.items():
                outname = self.prefix + '_' + key + '.nat'
                out = open(outname, 'w')
                out.write('ID\t' + self.prefix + '_' + key + '\n')
                for i in range(len(value)):
                    out.write(self.ids[i] + '\t' + value[i] + '\n')
                out.close()

        ## Main Write ##

        if f == 'attributes':
            _write_attributes()
        else:
            filename = file
            outfile = open(filename, "w")
            f = f.lower()

            if f == "nexus":
                _write_nexus(outfile, sep)
            elif f =="raxml":
                _write_raxml(outfile, sep)
            elif f == "plink":
                _write_plink(outfile, file, sep)
            elif f == "matrix":
                np.savetxt(filename, self.matrix, fmt='%.9f', delimiter=sep)
            elif f == "meta":
                _write_metadata(outfile, sep)
            elif f == "snp":
                _write_snpdata(outfile, sep)
            else:
                raise IOError("File format not supported.")

            outfile.close()
    
    def __str__(self):
        
        return ('-----------\nNumber of Individuals: %i\nNumber of SNPs: %i\nPloidy: %s\n-----------\n') % \
               (self.n, self.nSNP, self.ploidy)

#### Analysis Module ####

class Analysis:
    
    def __init__(self, data):
        
        self.data = data
        
    def getDistance(self, target='snps', distance='hamming'):

        print(get_time() + "\t" + 'Distance = ' + distance.upper())

        if self.data.filetype == 'dist':
            target = 'matrix'

        if target == 'matrix':
            matrix = np.array(self.data.matrix)
        else:
            # Convert alleles to numbers (e.g. A -> 1, B -> 2) for use in scipy.spatial.distance.pdist()
            allele_codes = {}
            for i in range(len(self.data.alleles)):
                allele_codes[self.data.alleles[i]] = int(i+1)
            allele_codes[self.data.missing] = 0   # missing can not be 1 to i

            snps = self.data.snps
            for a, code in allele_codes.items():
                snps[snps == a] = code

            matrix = snps

        if distance == 'asd':
            self.runPLINK(asd=True)
            self.data.readData(file=self.data.prefix + '_plink.mdist', f='matrix', sep=' ')
        else:
            matrix = sd.squareform(sd.pdist(matrix, distance))
            self.data.matrix = matrix

        self.data.matrices[distance] = self.data.matrix

        return matrix

    def runPLINK(self, qc_parameters={}, commandstring='', asd=False, quality=False):

        if self.data.ploidy == 'haploid':
            raise AttributeError('Haploid genotypes not supported for PLINK.')

        if commandstring:
            subprocess.call(commandstring)
        else:
            self.data.writeData(file=self.data.prefix + '_plink_in.ped', f='plink')

            if quality and qc_parameters:
                command = ['plink', '--noweb', '--file', self.data.prefix + '_plink_in']
                for key, value in qc_parameters.items():
                    command.append(key)
                    command.append(str(value))
                command.append('--recode')
                command.append('--out')
                command.append(self.data.prefix + '_plink_qc')

                subprocess.call(command, stdout=subprocess.DEVNULL)
                if os.path.exists(self.data.prefix + '_plink_qc.ped'):
                    self.data.readData(file=self.data.prefix + '_plink_qc.ped', f='plink', sep=' ')

            if asd:
                subprocess.call(['plink', '--noweb', '--file', self.data.prefix + '_plink_in', '--cluster', '--distance-matrix',
                                 '--out', self.data.prefix + '_plink'], stdout=subprocess.DEVNULL)

    def updateNodeAttributes(self, attribute_file):

        if os.path.isfile(self.data.prefix + '_plink_qc.irem'):
            infile = open(self.data.prefix + '_plink_qc.irem')
            to_remove = [line.strip().split()[1] for line in infile]
            infile.close()

            infile = open(attribute_file)
            outname = attribute_file.split('.')[0] + '_qc.csv'
            outfile = open(outname, 'w')
            for line in infile:
                content = line.strip().split(',')
                if content[0] not in to_remove:
                    outfile.write(line)
            infile.close()
            outfile.close()

            self.data.readData(file=outname, f='attributes', sep=',')

    def runNetView(self, tree=True, start=10, stop=40, step=10, algorithm='auto'):

        print(get_time() + "\t" + "Minimum Spanning Tree = " + str(tree).upper())
        print(get_time() + "\t" + "Nearest Neighbour = " + algorithm.upper())
        print(get_time() + "\t" + "k = " + str(start) + " - " + str(stop) + ' (by ' + str(step) + ')')
        print(get_time() + "\t" + "---------------------------------")

        self.data.netview_runs += 1
        
        matrix = self.data.matrix

        if tree:
            mst = csg.minimum_spanning_tree(matrix)
            mst = mst.toarray()
            self.data.networks['mst_' + str(self.data.netview_runs)] = mst
            mst = mst + mst.T
        else:
            mst = None

        pool = mp.Pool()
        networks = [pool.apply_async(netview, args=(matrix, k, mst, algorithm, tree,), callback=netview_callback)
                    for k in range(start, stop+1, step)]
        pool.close()
        pool.join()

        for item in networks:
            result = item.get()
            edges_array = result[1]
            edges = result[1].tolist()
            mst_edges = result[4].tolist()

            self.data.networks['netview_k' + str(result[0]) + '_' + str(self.data.netview_runs)] = result[1:]

            filename = self.data.prefix + '_netview_k'  + str(result[0]) +\
                       "_" + str(self.data.netview_runs) + '.edges'

            out = open(filename, "w")
            out.write('Source\tTarget\tDistance\tMST\n')
            for i in range(len(edges)):
                out.write(str(self.data.ids[edges[i][0]]) + "\t" + str(self.data.ids[edges[i][1]]) +
                          "\t" + str(matrix[edges[i][0], edges[i][1]]))
                if tree:
                    if edges[i] in mst_edges:
                        out.write('\t' + 'red\n')
                    else:
                        out.write('\t' + 'grey\n')
                else:
                    out.write("\n")

            if not tree:
                singletons = np.setdiff1d(np.arange(self.data.n), edges_array.flatten()).tolist()
                if singletons:
                    for node in singletons:
                        out.write(str(node) + '\n')
            out.close()

if __name__ == '__main__':
    main()
