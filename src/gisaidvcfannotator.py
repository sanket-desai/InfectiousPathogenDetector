#!/usr/bin/python3
'''

Author	  : Sonal Rashmi
Date		: 17/07/2020
Description : Clade Assessment 
'''

import pysam
import re
import sys
import os
import numpy as np
from gisaidvcfparser import *
from gisaidmetadataparser import *
from scipy.spatial import distance
from globals import *
import pandas as pd


class CustomError(Exception):
	def __int__(self,value):
		self.value=value
	def __str__(self):
		return repr(self.value)
class GisaidVcfAnnotator(object):
	def __init__(self, VcfFileDict=None):
		if VcfFileDict:
			self.GisaidVcfParserObj=GisaidVcfParser(GlobalVar.vcfdatabasefile_) 
			self.GisaidMetadataParserObj=GisaidMetadataParser(GlobalVar.gisaidmetadatafile_)
			self.key_list=[]
			self.annotation_dict_list=[]
			self.sample_list=[]
			for sample in VcfFileDict:
				self.sample_list.append(sample)
				VcfFile=VcfFileDict[sample]
				if os.path.exists(VcfFile) and os.stat(VcfFile).st_size > 0:
					self.vcf = pysam.VariantFile(VcfFile)
					self.key_temp=[]
					self.annotation_dict={}
					for rec in self.vcf.fetch():
						if rec.chrom.split('.')[0] == "NC_045512":
							key=str(str(rec.chrom.split('.')[0])+"_"+str(rec.pos)+"_"+str(rec.ref)+"_"+str(','.join(str(x) for x in rec.alts)))
							self.annotation_dict[key]=[rec.info['ANN'][0].split('|')[1],rec.info['ANN'][0].split('|')[3],rec.info['ANN'][0].split('|')[4],rec.info['ANN'][0].split('|')[10]]	
							self.key_temp.append(key)
					self.annotation_dict_list.append(self.annotation_dict)
					self.key_list.append(self.key_temp)
				else:
					self.annotation_dict_list.append(self.annotation_dict)
					self.key_list.append(self.key_temp)
		else:
			sys.exit(0)
	
	def total_variants_list(self):
		variants_count_list=[]
		for vcf_sample_key in self.key_list:
			variants_count_list.append(len(vcf_sample_key))
		return variants_count_list 
	
	def get_novel_variant_key_list(self):
		novel_variant_key_list=[]
		for vcf_sample_key in self.key_list:
			novel_key_temp=[]
			for key in vcf_sample_key:
				if not self.GisaidVcfParserObj.is_genomic_key_present(key):
					novel_key_temp.append(key)
			novel_variant_key_list.append(novel_key_temp)
		return novel_variant_key_list
			
	def get_novel_variant_annotation_dict_list(self):
		novel_variant_dict_list=[]
		novel_variant_key_list=self.get_novel_variant_key_list()
		cnt=0
		for vcf_sample_key in self.key_list:
			novel_variant_dict={}
			novel_variant_key=novel_variant_key_list[cnt]
			if len(novel_variant_key) > 0:
				for key in novel_variant_key:
					temp=self.annotation_dict_list[cnt][key]
					novel_variant_dict[key]=temp					
			cnt=cnt+1
			novel_variant_dict_list.append(novel_variant_dict)
		return novel_variant_dict_list	
			
	def get_sample_mut_profile(self):
		sample_mut_profile_list=[]
		for vcf_sample_key in self.key_list:
			sample_mut_profile = [0] * len(self.GisaidVcfParserObj.get_genomic_key())
			for key in vcf_sample_key:
				if self.GisaidVcfParserObj.is_genomic_key_present(key):
					key_index=self.GisaidVcfParserObj.get_genomic_key().index(key)		
					sample_mut_profile[key_index]=1
			sample_mut_profile_list.append(sample_mut_profile)
		return sample_mut_profile_list
	
	def get_gisaid_for_all_samples(self):
		gisaid_list_all_sample_list=[]
		sample_mut_profile_list=self.get_sample_mut_profile()
		cnt=0
		for vcf_sample_key in self.key_list:
			euclidean_gisaid_list=[]
			gisaid_list=[]
			gisaid_id_list=[]
			gisaid_mut_profile = [] 
			for key in vcf_sample_key:
				if self.GisaidVcfParserObj.is_genomic_key_present(key):
					gsaid_complete_list=self.GisaidVcfParserObj.get_gisaid_given_genomic_key(key)
					for id in gsaid_complete_list:
						gisaid_id_list.append(id)
						mut_profile=self.GisaidVcfParserObj.get_mutational_profile_given_gisaid(id)
						gisaid_mut_profile.append(mut_profile)
			sample_mut_profile=sample_mut_profile_list[cnt]
			if len(gisaid_id_list) > 0:
				euclidean_distance_list=[]
				for mut in gisaid_mut_profile:
					euclidean_distance_list.append(distance.euclidean(sample_mut_profile,mut))
				min_euclidean_value=min(euclidean_distance_list)
				id_index=list(np.where(euclidean_distance_list==min_euclidean_value)[0])
				for index in id_index:
					id=gisaid_id_list[index]
					gisaid_list.append(id)
				euclidean_gisaid_list=[min_euclidean_value,gisaid_list]	
			gisaid_list_all_sample_list.append(euclidean_gisaid_list)
			cnt=cnt+1
		return gisaid_list_all_sample_list
		
	def get_gisaid_clade_for_sample(self):
		gisaid_clade_list=[]
		sample_list_gisaid_id=self.get_gisaid_for_all_samples()
		cnt=0
		for gisaid_sample_list in sample_list_gisaid_id:
			id_clade_ann_list= [None] * 4
			if len(gisaid_sample_list) > 0:
				distance=gisaid_sample_list[0]
				gisaid_list=gisaid_sample_list[1]
				clade_dict={}
				if len(gisaid_list) > 0: 
					for gisaid_id in gisaid_list:				
						gisaid_clade=self.GisaidMetadataParserObj.get_gisaid_clade_for_gisaid_id(gisaid_id)
						pangolin_lineage=self.GisaidMetadataParserObj.get_pangolin_lineage_for_gisaid_id(gisaid_id)
						key=str(gisaid_clade)+","+str(pangolin_lineage)
						if key not in clade_dict:
							clade_dict[key]=[]
							clade_dict[key].append(gisaid_id)
						else:
							clade_dict[key].append(gisaid_id)
					clade_list=[]
					lineage_list=[]
					id_list=[]
					
					for clade in clade_dict:
						clade_list.append(clade.split(',')[0])
						lineage_list.append(clade.split(',')[1])
						id_list.append(clade_dict[clade][0]+" - "+str(len(clade_dict[clade])))
					id_clade_ann_list=[str(';'.join(str(x) for x in id_list)),str(';'.join(str(x) for x in clade_list)),str(';'.join(str(x) for x in lineage_list)),distance]		
			gisaid_clade_list.append(id_clade_ann_list)
		return gisaid_clade_list
	
	def get_pandas_df_novel_variant(self):
		df=pd.DataFrame()
		col=["Sample","Genome","Position","Reference","Altered","Consequence","Gene","Transcript","Protein Change"]
		novel_mut_anno_dict_list=self.get_novel_variant_annotation_dict_list()
		cnt=0
		data=[]
		for sample in self.sample_list:
			if len(novel_mut_anno_dict_list[cnt]) > 0:
				for key in novel_mut_anno_dict_list[cnt]:
					record=[str(sample),key.split("_")[0]+"_"+key.split("_")[1],key.split("_")[2],key.split("_")[3],key.split("_")[4],novel_mut_anno_dict_list[cnt][key][0],novel_mut_anno_dict_list[cnt][key][1],novel_mut_anno_dict_list[cnt][key][2],novel_mut_anno_dict_list[cnt][key][3]]
					data.append(record)
			else:
				record=[str(sample),"","","","","","","",""]
				data.append(record)
			cnt=cnt+1			
		df = pd.DataFrame(data, columns = col)
		return df
		
	def get_pandas_df_clade_assessment(self):
		df=pd.DataFrame()
		col=["Sample","Number Of Variants","Related Strain","GISAID Clade","Pangolin Lineage","Euclidean Distance"]
		total_variant_list=self.total_variants_list()
		gisaid_clade_list=self.get_gisaid_clade_for_sample()
		cnt=0
		data=[]
		for sample in self.sample_list:
			if gisaid_clade_list[cnt][0]:
				record=[str(sample),str(total_variant_list[cnt]),str(gisaid_clade_list[cnt][0]),str(gisaid_clade_list[cnt][1]),str(gisaid_clade_list[cnt][2]),str(gisaid_clade_list[cnt][3])]
			else:
				record=[str(sample),"","","","",""]
			data.append(record)
			cnt=cnt+1
		df = pd.DataFrame(data, columns = col)
		return df
			
