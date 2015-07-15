#!/usr/bin/env python
# NetView P v.0.7.1
# Dependencies: PLINK
# Eike Steinig
# Zenger Lab, JCU
# https://github.com/esteinig/netview

import os
import time
import json
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
    dat.project = commands.arg_dict['project']
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
    elif commands.arg_dict['nexus']:
        dat.filetype = 'nexus'
        dat.readData(commands.arg_dict['data_file'], f='nexus', sep=commands.arg_dict['sep'])
    elif commands.arg_dict['snps']:
        dat.filetype = 'raxml'
        dat.readData(commands.arg_dict['data_file'], f='raxml', sep=commands.arg_dict['sep'])
    else:
        dat.filetype = 'dist'
        dat.readData(commands.arg_dict['data_file'], f='matrix', sep=commands.arg_dict['sep'])

    dat.readData(commands.arg_dict['attribute_file'], f='attributes', sep=',')

    for stored in dat.meta_data.values():
        if len(stored) != dat.n:
            print('\nError. N in Data != N in Attribute File.')
            exit(1)

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

    if commands.arg_dict['qc'] and pipeline.data.filetype != 'dist':

        qc_params = {'--mind': commands.arg_dict['mind'],
                     '--geno': commands.arg_dict['geno'],
                     '--maf': commands.arg_dict['maf'],
                     '--hwe': commands.arg_dict['hwe']}

        pipeline.runPLINK(qc_parameters=qc_params, quality=True)
        pipeline.updateNodeAttributes(commands.arg_dict['attribute_file'])

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
                            algorithm=commands.arg_dict['algorithm'], edges=commands.arg_dict['edges'],
                            html=commands.arg_dict['web'])

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
    d3_path = os.path.realpath(os.path.join(project_path, 'd3'))

    if os.path.exists(project_path):
        shutil.rmtree(project_path)

    architecture = [project_path, plink_path, network_path, other_path, node_path, d3_path]

    for directory in architecture:
        try:
            os.makedirs(directory)
        except OSError:
            if not os.path.isdir(directory):
                raise

    for name in os.listdir(cwd):
        pathname = os.path.join(cwd, name)
        if os.path.isfile(pathname):
            if name.endswith('.edges'):
                shutil.move(pathname, network_path)
            elif name.endswith('.dist'):
                shutil.move(pathname, other_path)
            elif name.endswith('.nat'):
                shutil.move(pathname, node_path)
            elif name.startswith(prefix + '_plink'):
                shutil.move(pathname, plink_path)
            elif name.endswith('_qc.csv'):
                shutil.move(pathname, other_path)
            elif name.endswith('.json') or name.endswith('.html'):
                shutil.move(pathname, d3_path)
            elif name.startswith(prefix + '_plink_in'):
                os.remove(pathname)

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

    if tree:
        mst_edges = np.argwhere(adjacency < 1)
    else:
        mst_edges = np.array([])

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
        data_type.add_argument('-n', dest='nexus', action='store_true',
                               help="Nexus format from SPANDx")
        data_type.add_argument('-r', dest='raxml', action='store_true',
                               help="RAxML format from SPANDx")
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

        self.parser.add_argument('--html', dest='web', action='store_true', default=True,
                                 help="Generate D3/JSON graphs (ON)")
        self.parser.add_argument('--edges', dest='edges', action='store_true', default=False,
                                 help="Generate graphs as edge files (OFF)")

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

        self.project = "project"
        self.prefix = "prefix"
        self.ploidy = 'diploid'
        self.missing = "0"
        
        self.n = 0
        self.nSNP = 0
        self.ids = []  # IDs

        self.alleles = []
        self.snps = np.array([])  # Array (N x SNPs)
        self.biodata = []  # List/Alignment of BioPython SeqRecords
        
        self.meta_data = {}
        self.snp_data = {}
        
        self.matrices = {}
        self.networks = {}

        self.matrix = np.array([])  # Current Matrix
        
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
                if matrix:
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

        def _read_attributes(file):
            content = [line.strip().split(',') for line in file]
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
            _read_attributes(infile)
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

        def _write_graph_json():

            col_dict = {'dimgray': '#696969', 'olive': '#808000', 'burlywood': '#deb887', 'darkgreen': '#006400',
                    'navy': '#000080', 'white': '#ffffff', 'violet': '#ee82ee', 'darkblue': '#00008b',
                    'steelblue': '#4682b4', 'deepskyblue': '#00bfff', 'tan': '#d2b48c', 'rebeccapurple': '#663399',
                    'honeydew': '#f0fff0', 'slategray': '#708090', 'powderblue': '#b0e0e6', 'palevioletred': '#db7093',
                    'chocolate': '#d2691e', 'coral': '#ff7f50', 'azure': '#f0ffff', 'peru': '#cd853f',
                    'springgreen': '#00ff7f', 'darkorange': '#ff8c00', 'mediumvioletred': '#c71585',
                    'mediumaquamarine': '#66cdaa', 'darkmagenta': '#8b008b', 'mediumslateblue': '#7b68ee',
                    'mediumseagreen': '#3cb371', 'crimson': '#dc143c', 'gainsboro': '#dcdcdc', 'darkgray': '#a9a9a9',
                    'plum': '#dda0dd', 'forestgreen': '#228b22', 'seagreen': '#2e8b57', 'teal': '#008080',
                    'gold': '#ffd700', 'dodgerblue': '#1e90ff', 'lightpink': '#ffb6c1', 'papayawhip': '#ffefd5',
                    'orchid': '#da70d6', 'black': '#000000', 'cornflowerblue': '#6495ed', 'lightyellow': '#ffffe0',
                    'goldenrod': '#daa520', 'purple': '#800080', 'khaki': '#f0e68c', 'aquamarine': '#7fffd4',
                    'lightskyblue': '#87cefa', 'fuchsia': '#ff00ff', 'mediumblue': '#0000cd', 'sandybrown': '#f4a460',
                    'moccasin': '#ffe4b5', 'darkslategray': '#2f4f4f', 'cornsilk': '#fff8dc', 'lightcyan': '#e0ffff',
                    'darkolivegreen': '#556b2f', 'silver': '#c0c0c0', 'lightgoldenrodyellow': '#fafad2',
                    'navajowhite': '#ffdead', 'turquoise': '#40e0d0', 'rosybrown': '#bc8f8f', 'antiquewhite': '#faebd7',
                    'thistle': '#d8bfd8', 'lightcoral': '#f08080', 'floralwhite': '#fffaf0', 'indianred': '#cd5c5c',
                    'ghostwhite': '#f8f8ff', 'blue': '#0000ff', 'snow': '#fffafa', 'orangered': '#ff4500',
                    'darkred': '#8b0000', 'greenyellow': '#adff2f', 'ivory': '#fffff0', 'mediumorchid': '#ba55d3',
                    'lawngreen': '#7cfc00', 'lightsalmon': '#ffa07a', 'lightgray': '#d3d3d3',
                    'lightslategray': '#778899', 'mediumpurple': '#9370db', 'darkcyan': '#008b8b', 'tomato': '#ff6347',
                    'lightsteelblue': '#b0c4de', 'darkseagreen': '#8fbc8f', 'aqua': '#00ffff', 'olivedrab': '#6b8e23',
                    'darkgoldenrod': '#b8860b', 'darkorchid': '#9932cc', 'seashell': '#fff5ee', 'skyblue': '#87ceeb',
                    'blanchedalmond': '#ffebcd', 'beige': '#f5f5dc', 'darkturquoise': '#00ced1', 'slateblue': '#6a5acd',
                    'red': '#ff0000', 'lavender': '#e6e6fa', 'hotpink': '#ff69b4', 'yellowgreen': '#9acd32',
                    'cyan': '#00ffff', 'firebrick': '#b22222', 'lemonchiffon': '#fffacd', 'darksalmon': '#e9967a',
                    'sienna': '#a0522d', 'mediumturquoise': '#48d1cc', 'salmon': '#fa8072', 'green': '#008000',
                    'lightgreen': '#90ee90', 'deeppink': '#ff1493', 'palegoldenrod': '#eee8aa', 'orange': '#ffa500',
                    'wheat': '#f5deb3', 'lime': '#00ff00', 'lavenderblush': '#fff0f5', 'brown': '#a52a2a',
                    'blueviolet': '#8a2be2', 'magenta': '#ff00ff', 'lightseagreen': '#20b2aa', 'mistyrose': '#ffe4e1',
                    'saddlebrown': '#8b4513', 'midnightblue': '#191970', 'mediumspringgreen': '#00fa9a',
                    'cadetblue': '#5f9ea0', 'paleturquoise': '#afeeee', 'palegreen': '#98fb98', 'pink': '#ffc0cb',
                    'darkkhaki': '#bdb76b', 'oldlace': '#fdf5e6', 'whitesmoke': '#f5f5f5', 'royalblue': '#4169e1',
                    'gray': '#808080', 'lightblue': '#add8e6', 'maroon': '#800000', 'peachpuff': '#ffdab9',
                    'darkslateblue': '#483d8b', 'linen': '#faf0e6', 'limegreen': '#32cd32',
                    'mintcream': '#f5fffa', 'chartreuse': '#7fff00', 'yellow': '#ffff00', 'indigo': '#4b0082',
                    'bisque': '#ffe4c4', 'aliceblue': '#f0f8ff', 'darkviolet': '#9400d3'}

            if self.networks.keys() == '':
                print('No networks to write to JSON.')

            templ_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')
            templ_file = open(os.path.join(templ_path, 'fd_network.html'))
            templ_str = templ_file.read()
            templ_file.close()

            for network, properties in self.networks.items():

                json_name = self.prefix + '_' + network + '.json'
                html_name = self.prefix + '_' + network + '.html'
                edges = properties[1]
                weights = properties[2]

                node_array = []
                for i in range(len(self.ids)):
                    if 'lat' in self.meta_data.keys() and 'lon' in self.meta_data.keys():
                        node_array.append({'name': self.ids[i], 'group': self.meta_data['pop'][i], 'color':
                                           col_dict[self.meta_data['col'][i]], 'lon': self.meta_data['lon'][i],
                                           'lat': self.meta_data['lat'][i]})
                    else:
                        node_array.append({'name': self.ids[i], 'group': self.meta_data['pop'][i], 'color':
                                            col_dict[self.meta_data['col'][i]]})

                edge_array = []
                for i in range(len(edges)):
                    if 'lat' in self.meta_data.keys() and 'lon' in self.meta_data.keys():
                        edge_array.append({'source': int(edges[i, 0]), 'target': int(edges[i, 1]), 'value':
                                    float(weights[i]), 'slon': self.meta_data['lon'][edges[i, 0]], 'slat':
                                           self.meta_data['lat'][edges[i, 0]], 'tlon': self.meta_data['lon'][edges[i, 1]],
                                           'tlat': self.meta_data['lat'][edges[i, 1]]})
                    else:
                        edge_array.append({'source': int(edges[i, 0]), 'target': int(edges[i, 1]), 'value':
                                    float(weights[i])})

                json_file = open(json_name, 'w')
                json_file.write(json.dumps({'nodes': node_array, 'links': edge_array, }, sort_keys=True, indent=2))
                json_file.close()

                if self.ploidy == 'diploid':
                    nsnps = self.nSNP//2
                else:
                    nsnps = self.nSNP

                html_file = open(html_name, 'w')
                html = templ_str.replace('template.json', '"' + json_name + '"')
                html = html.replace('temp_n', str(self.n))
                html = html.replace('temp_snp', str(nsnps))
                html = html.replace('temp_k', str(properties[0]))
                html = html.replace('temp_project', str(self.project))
                html_file.write(html)
                html_file.close()

        def _write_graph_edges():

            if self.networks.keys() == '':
                print(get_time() + '\t' + 'Warning: No networks to write to JSON.')

            for network, properties in self.networks.items():
                filename = network + '.edges'

                edges = properties[1].tolist()
                weights = properties[2].tolist()
                mst_edges = properties[4].tolist()

                out = open(filename, 'w')
                out.write('Source\tTarget\tDistance\tMST\n')

                for i in range(len(edges)):
                    out.write(str(self.ids[edges[i][0]]) + "\t" + str(self.ids[edges[i][1]]) +
                          "\t" + str(weights[i]))

                    if len(mst_edges) > 0:
                        if edges[i] in mst_edges:
                            out.write('\t' + 'red\n')
                        else:
                            out.write('\t' + 'grey\n')
                    else:
                        out.write("\n")

                if len(mst_edges) == 0:
                    singletons = np.setdiff1d(np.arange(self.n), edges.flatten()).tolist()
                    if singletons:
                        for node in singletons:
                            out.write(str(self.ids[node]) + '\n')
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
            elif f == "json":
                _write_graph_json()
            elif f == 'edges':
                _write_graph_edges()
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

    def runNetView(self, tree=True, start=10, stop=40, step=10, algorithm='auto', edges=False, html=True):

        print(get_time() + "\t" + "Minimum Spanning Tree = " + str(tree).upper())
        print(get_time() + "\t" + "Nearest Neighbour = " + algorithm.upper())
        print(get_time() + "\t" + "k = " + str(start) + " - " + str(stop) + ' (by ' + str(step) + ')')
        print(get_time() + "\t" + "---------------------------------")

        self.data.netview_runs += 1
        
        matrix = self.data.matrix

        if tree:
            mst = csg.minimum_spanning_tree(matrix)
            mst = mst.toarray()
            #self.data.networks[self.data.prefix + 'mst_' + str(self.data.netview_runs)] = mst
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
            self.data.networks['netview_k' + str(result[0]) + '_' + str(self.data.netview_runs)] = result
        print(get_time() + "\t" + "---------------------------------")

        if html:
            print(get_time() + "\t" + "Out = JSON")
            self.data.writeData(f='json')
        if edges:
            self.data.writeData(f='edges')
            print(get_time() + "\t" + "Out = Edges")

main()
