#!/usr/bin/python3
'''

Author	  : Sonal Rashmi
Date		: 16/07/2020
Description : IPD Command Line Interface with two subparser long and short read. 

'''
import argparse
import os
import sys
from commandlineparser import *
from ipdshortread import *
from ipdlongread import *
import pathlib

#Directory and file check
def file_choices(choices,fname):
	ext = ''.join(pathlib.Path(fname).suffixes)
	if ext.lower() not in choices:
		print(ext)
		raise argparse.ArgumentTypeError('File must have proper extension')
	return fname
	
def dir_path(path):
	if os.path.isdir(path):
		return path
	else:
		raise NotADirectoryError(path)	

#define parser and subparsers  
def main():
	parent_parser = argparse.ArgumentParser(add_help=False)
	parent_parser.add_argument('-p',dest='prefix',type=str,required=True,help="set project name")
	parent_parser.add_argument('-t',dest='thread',type=int,default=4,required=False,help="set threads (default = 4)")
	parent_parser.add_argument('-m',dest='molecular_type',type=str,choices=['dna','rna'],required=True,help="set molecular type (DNA or RNA)")
	parent_parser.add_argument('-o',dest='outdir',type= lambda path:dir_path(path),required=True,default=os.getcwd(),help="set Output Directory")
	
	parser = argparse.ArgumentParser(add_help=False) 
	subparsers = parser.add_subparsers(dest='parser_name')

	# subcommand long parser																   
	parser_SingleEndSample = subparsers.add_parser('long', parents = [parent_parser])						  
	parser_SingleEndSample.add_argument('-i',dest='input_file',required=True,type=lambda s:file_choices((".fastq",".fastq.gz",".fq.gz",".fq"),s),help="Fastq file (With Complete Path)")					 

	# subcommand short parser																   
	parser_PairedEndSample = subparsers.add_parser('short', parents = [parent_parser])
	parser_PairedEndSample.add_argument('-i',dest='input_files',nargs='+',required=True,type=lambda s:file_choices((".fastq",".fastq.gz",".fq.gz",".fq"),s),help="Fastq file/files (With Complete Path)")
	
	args = parser.parse_args()													   
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	else:
		if args.parser_name == "long":
			runmode="SingleEndSample"
			sampleinfo=None
			fastq2 = None
			return_list=[runmode, args.outdir, args.prefix, args.molecular_type, args.thread, sampleinfo, args.input_file, fastq2]
			Parser=CommandLineInterfaceParser(return_list)
			maplist=Parser.map_single_end_sample()
			for inmap in maplist:				
				try:
					i=IPDLongRead(inmap)
				except TypeError:
					parser.print_help()
					#print(inmap)
					sys.exit(0)	
		elif args.parser_name == "short":
			input_file_len=len(args.input_files)
			if input_file_len == 1:
				runmode="SingleEndSample"
				sampleinfo=None
				fastq2 = None
				return_list=[runmode, args.outdir, args.prefix, args.molecular_type, args.thread, sampleinfo, args.input_files[0], fastq2]
				Parser=CommandLineInterfaceParser(return_list)
				maplist=Parser.map_single_end_sample()
				for inmap in maplist: 					
					try:
						i=IPDShortRead(inmap)
					except TypeError:
						parser.print_help()
						#print(inmap)
						sys.exit(0)
			elif input_file_len == 2:
				runmode="PairedEndSample"
				sampleinfo = None
				return_list=[runmode, args.outdir, args.prefix, args.molecular_type, args.thread, sampleinfo, args.input_files[0], args.input_files[1]]
				Parser=CommandLineInterfaceParser(return_list)
				maplist=Parser.map_paired_sample()
				for inmap in maplist: 					
					try:
						i=IPDShortRead(inmap)
					except TypeError:
						parser.print_help()
						#print(inmap)
						sys.exit(0)
			else:
				print("Invalid Input")
				sys.exit()
		else:
			print("Invalid Option")
			sys.exit()
if __name__ =="__main__":
	main()
