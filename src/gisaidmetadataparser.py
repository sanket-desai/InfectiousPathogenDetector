import pandas as pd
import re
import sys
import os

class CustomError(Exception):
	def __int__(self,value):
		self.value=value
	def __str__(self):
		return repr(self.value)
class GisaidMetadataParser(object):
	def __init__(self, MetadataFile=None):
		if os.path.exists(MetadataFile):
			self.GisaidMetadataPandasDf=pd.read_csv(MetadataFile,sep="\t",engine='python')
		else:
			sys.exit(0)
	
	def get_country_for_gisaid_id(self,id):
		country=None
		if id:
			if id in self.GisaidMetadataPandasDf.values:
				country=self.GisaidMetadataPandasDf.loc[self.GisaidMetadataPandasDf['gisaid_epi_isl'] == id]["country"].to_string(index=False).lstrip()
		return country

	def get_gisaid_clade_for_gisaid_id(self,id):
		gisaid_clade=None
		if id:
			if id in self.GisaidMetadataPandasDf.values:
				gisaid_clade=self.GisaidMetadataPandasDf.loc[self.GisaidMetadataPandasDf['gisaid_epi_isl'] == id]["GISAID_clade"].to_string(index=False).lstrip()
		return gisaid_clade

	def get_pangolin_lineage_for_gisaid_id(self,id):
		pangolin_lineage=None
		if id:
			if id in self.GisaidMetadataPandasDf.values:
				pangolin_lineage=self.GisaidMetadataPandasDf.loc[self.GisaidMetadataPandasDf['gisaid_epi_isl'] == id]["pangolin_lineage"].to_string(index=False).lstrip()
		return pangolin_lineage
