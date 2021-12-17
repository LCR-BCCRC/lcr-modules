#!/usr/bin/python
# -*- coding: utf-8 -*-
#The script builds an evolutionary history of clones by sharing mutations between samples.

import os
import re
import numpy 
import random 
import sys, getopt
from Bio import Phylo
from io import StringIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#24 chromosome
all_chr=[]
for i in range(1,24):
	all_chr.append(str(i))
all_chr.append('X')
all_chr.append('Y')


#parameter

opts,args = getopt.getopt(sys.argv[1:], "i:o:",["AF=", "drivers=", "input=", "output="])
def get_driver_gene_path():
    '''
    return the driver gene path
    '''
    driver_gene_path = os.path.dirname(os.path.realpath(__file__))+os.path.sep+'299_driverMutationList_Cell_2018.txt'
    return driver_gene_path


AF=''
file_path=''
file_path_out=''
for para,value in opts:
	if para in ('--AF'):
		AF=value
	if para in ('--drivers'): 
		drivers=value
	if para in ('-i','--input'):
		file_path=value
	if para in ('-o','--output'):
		file_path_out=value
print('AF: '+str(AF))
if file_path[-len(os.path.sep):]==os.path.sep:
	file_path=file_path[:-len(os.path.sep)]
if drivers: 
	print('drivers: '+drivers)
else: 
	drivers = get_driver_gene_path()
	print('drivers: '+drivers)
print('file_path: '+str(file_path))
if file_path_out[-len(os.path.sep):]==os.path.sep:
	file_path_out=file_path_out[:-len(os.path.sep)]
print('file_path_out: '+str(file_path_out))

def fact(n):
    if n == 0:
        return 1
    else:
        return n*fact(n-1)
def Cmn(n,m):
    return fact(n)/(fact(n-m)*fact(m))

#map driver gene:
# @ a is a list including mutation which is used to map driver gene
# @ b =='FALSE', then output driver gene without mutation
# @ b =='TRUE', then output driver gene with mutation
def search_gene(a,b):
	global count_gene
	count_gene=0
	global gene_list
	gene_list=[]
	for item in a:#read a mutation
		if item in all_mutation.keys():
			for gene in all_mutation[item]:
				if gene in driver_gene.keys():
					if b=='FALSE':#do not output driver gene with mutation
						gene=gene
					if b=='TRUE':
						gene=gene+':'+item
					if gene not in gene_list:
						gene_list.append(gene)
	count_gene=len(gene_list)
	gene_list.append(count_gene)
	return gene_list							
									
#.................................................find best branch..........................
def xunhuan(a):
	num=1
	global first_1
	first_1={}
	global first_3
	first_3={}
	global first_4
	first_4={}
	first_6={}
	while num<=len(a)/2:
		struct_pattern=str(num)+'-'+str(len(a)-num)
		DIC={}#including structure info: Cmn information for a specific branch split
		cishu=int(Cmn(len(a),num))
		while len(DIC.keys())<cishu:
			sel=random.sample(a,num)
			sel.sort()
			if num>1:
				da=','
				for Item in sel:
					da=da+','+Item
				if da[2:] not in DIC.keys():
					DIC.update({da[2:]:[]})
			if num==1:
				DIC.update({sel[0]:[]})
		first_2={}#struct info
		for item in DIC.keys():
			a_1=[]#first =a, then remove item from [a]. 
			jiaoji_2=[]
			jiaoji_1=[]
			for sample in a:
				a_1.append(sample)
			split=item.split(',')
			split.sort()
			for item5 in split:
				a_1.remove(item5)
			if len(split)==1:
				jiaoji_1=dic[item]#sample mutation which is removed
			else:
				NUM=0
				for item6 in split:
					NUM=NUM+1
					if NUM==1:
						jiaoji_1=dic[item6]
					jiaoji_1=[val for val in jiaoji_1 if val in dic[item6]]#sample mutation which are removed

			i=0
			no_item=','
			a_1.sort()
			for item1 in a_1:
				no_item=no_item+','+item1
				i=i+1
				if i==1:
					jiaoji_2=dic[item1]
				jiaoji_2 = [val for val in jiaoji_2 if val in dic[item1]]#sample mutation which are retained; also is big group

			struct_name=no_item[2:]
			struct_share_mut=len(jiaoji_2)
			num_1=str(len(jiaoji_1))+str(search_gene(jiaoji_1,'FALSE'))
			first_3.update({item:num_1})
			first_4.update({item:jiaoji_1})
			num_2=str(len(jiaoji_2))+str(search_gene(jiaoji_2,'FALSE'))
			first_3.update({struct_name:num_2})
			first_4.update({struct_name:jiaoji_2})
			#if struct: n/2,n/2: we should retained group info which mutation is small
			if num==len(a)/2:
				if len(jiaoji_1)>len(jiaoji_2):
					struct_name=no_item[2:]
					struct_share_mut=len(jiaoji_2)
				if len(jiaoji_1)<len(jiaoji_2):
					struct_name=item
					struct_share_mut=len(jiaoji_1)
				if len(jiaoji_1)==len(jiaoji_2):
					struct_name=item
					struct_share_mut=len(jiaoji_1)
			first_2.update({struct_name:struct_share_mut})
		sort_d=sorted(first_2.items(),key = lambda d:d[1],reverse=True) 
		count=0
		#print('max value and second value:')
		max_data=[]
		if num!=len(a)/2:
			for key,value in sort_d:
				count+=1
				if count==1:
					max_data_name=key
					max_data_value=value
					struct_pattern1=struct_pattern+':'+key
					first_6.update({struct_pattern1:value})
				if count<=2:
					max_data.append(value)
					#print(key,value)
				if count>1 and key!=max_data_name and value==max_data_value:
					struct_pattern1=struct_pattern+':'+key
					if struct_pattern1 not in first_6.keys():
						first_6.update({struct_pattern1:value})
		else:
			for key,value in sort_d:
				count+=1
				if count==1:
					max_data_name=key
					max_data_value=value
					max_data.append(value)
					struct_pattern1=struct_pattern+':'+max_data_name
					if struct_pattern1 not in first_6.keys():
						first_6.update({struct_pattern1:max_data_value})
					#print(key,value)
					chongfu={}#judge random chongfu 
					while len(chongfu.keys())<(fact(len(key.split(',')))-1):
						list1=key.split(',')
						random.shuffle(list1)
						#print(list1)
						fan_da=','
						for Item1 in list1:
							fan_da=fan_da+','+Item1
						if fan_da[2:]!=key:
							chongfu.update({fan_da[2:]:[]})
				if key != max_data_name and key not in chongfu.keys() and len(max_data)==1:
					max_data.append(value)
					#print(key,value)
				if count>1 and key not in chongfu.keys() and key!=max_data_name and value==max_data_value:
					chongfu2={}
					while len(chongfu2.keys())<(fact(len(key.split(',')))-1):
						list1=key.split(',')
						random.shuffle(list1)
						fan_da=','
						for Item1 in list1:
							fan_da=fan_da+','+Item1
						if fan_da[2:]!=key:
							chongfu2.update({fan_da[2:]:[]})
					chongfu_count=0
					for item_name in chongfu2.keys():
						struct_pattern1=struct_pattern+':'+item_name
						if struct_pattern not in first_6.keys():
							chongfu_count=chongfu_count+1
					if chongfu_count==len(chongfu.keys()):
						struct_pattern1=struct_pattern+':'+key
						if struct_pattern1 not in first_6.keys():
							first_6.update({struct_pattern1:value})
		if int(max_data[1])!=0:
			ratio=int(max_data[0])/int(max_data[1])
		else:
			ratio=int(max_data[0])
		struct_pattern2=struct_pattern+':'+max_data_name
		first_1.update({struct_pattern2:ratio})
		num=num+1
	sort_data=sorted(first_1.items(),key = lambda d:d[1],reverse=True) 
	#print(sort_data)
	global return_data
	return_data=[]
	count_data=0
	for key,value in sort_data:
		count_data=count_data+1
		if count_data==1:
			max_data_name=key
			max_data_value=value
		if value==max_data_value:
			return_data.append(key)
	for item in return_data:
		identify=item.split(':')[0]
		for item1 in first_6.keys():
			if identify in item1 :
				if item1 not in return_data:
					return_data.append(item1)
	if max_data_value==0:#ratio==0
		value1=[]
		return_data=[]
		num=divmod(len(a),2)[0]
		while num>0:
			for key in first_4.keys():
				if len(key.split(','))==num and len(first_4[key])!=0:
					value1.append(len(first_4[key]))
			if len(value1)!=0:
				break
			else:
				num=num-1
		max_data=value1[0]
		struct_pattern=str(num)+'-'+str(len(a)-num)
		for data in value1[1:]:
			if int(data)>max_data:
				max_data=int(data)
		for key in first_4.keys():
			if len(key.split(','))==num and len(first_4[key])==max_data:
				struct_pattern1=struct_pattern+':'+key
				return_data.append(struct_pattern1)
	return_data.sort()
	return(return_data)
#split tree
def split_tree(tree_label_split):
	new_tree_label=[]
	num=0
	count1=0
	count2=0
	index1=[]
	new_tree_label=[]
	while num <len(tree_label_split):
		if tree_label_split[num]=='(':
			count1=count1+1
		if tree_label_split[num]==')':
			count2=count2+1
		if count1==count2:
			if tree_label_split[num]==']':
				index1.append(num+1)
				count1=0
				count2=0
		num=num+1
	index_1=0
	for item in index1:
		new_tree_label.append(tree_label_split[index_1:item])
		index_1=item+1
	#print(new_tree_label)
	return(new_tree_label)	
	
#...........................................299 driver gene........................
driver_gene_path=drivers

print('driver_gene_path: '+driver_gene_path)
if driver_gene_path!='.':
	file=open(driver_gene_path,'r')
	global driver_gene
	driver_gene={}
	lines=file.readlines()
	for line in lines[1:]:
		driver_gene.update({line.rstrip():[]})  #uniq driver gene
else:
	print('can not read driver gene file')

#...................................................read file...and mutation................................

filelist=os.listdir(file_path)
print(filelist)
#find all patient
all_patient_ID=[]
for filename in filelist:
	patient_ID=filename.split('.')[0].split('_')[0]
	if patient_ID not in all_patient_ID:
		all_patient_ID.append(patient_ID)
global all_mutation
for patient in all_patient_ID:
	if AF !='':#one sample per file
		print('one sample per file')
		all_mutation={}
		print('patient: '+patient)
		dic={}#[sample_id: all mutation_ID]
		all_data=[]#all_sample_ID
		sample_shu=0#the number of all samples of this patient
		for filename in filelist:
			if filename.split('.')[0].split('_')[0]==patient:
				sample_shu=sample_shu+1
				sample_id=filename.split('.')[0].split('_')[1]
				if sample_id not in all_data:
					all_data.append(sample_id)
				else:
					print('sample_id is not uniq!')
				if sample_id not in dic.keys():
					dic.update({sample_id:[]})
				else:
					print('sample_ID is not uniq!')
				file=open(file_path+os.path.sep+filename,'r')
				lines=file.readlines()
				for line in lines[1:]:
					data=line.rstrip().split('\t')
					if int(data[1])+int(data[2])!=0:
						AF_data=int(data[1])/(int(data[1])+int(data[2]))
						if AF_data>=float(AF) and data[0].split('_')[0] in all_chr:
							if data[0] not in dic[sample_id] :
								dic[sample_id].append(data[0])
								if data[0] not in all_mutation.keys():
									all_mutation.update({data[0]:[]})
								for refer_gene in data[3].split(';'):
									if refer_gene not in all_mutation[data[0]]:
										all_mutation[data[0]].append(refer_gene)
							else:
								print(filename+' contains two or more identical mutation!')
					else:
						print(data[0]+': this sample have mutations that are uncovered!')
	else:
		# 0-1 matrix
		print('0-1 matrix')
		all_mutation={}
		print('patient: '+patient)
		dic={}#[sample_id: all mutation_ID]
		all_data=[]#all_sample_ID
		sample_shu=0#the number of all samples of this patient
		sample_index={}#[index:smaple_name]
		file=open(file_path+os.path.sep+patient+'.txt','r')
		lines=file.readlines()
		data_0=lines[0].rstrip().split('\t')
		for item in data_0[1:-1]:#sample_id
			sample_shu=sample_shu+1
			if item not in dic.keys() :
				dic.update({item:[]})
				all_data.append(item)
				sample_index.update({sample_shu:item})
			else:
				print('sample_id is not uniq!')
		#print(sample_index)
		#print(dic)
		for line in lines[1:]:
			data=line.rstrip().split('\t')
			for i in range(1,len(data)-1):
				if int(data[i])==1 and data[0] not in dic[sample_index[i]]:
					dic[sample_index[i]].append(data[0])
			if data[0] not in all_mutation.keys():
				all_mutation.update({data[0]:[]})
				for refer_gene in data[-1].split(';'):
					if refer_gene not in all_mutation[data[0]]:
						all_mutation[data[0]].append(refer_gene)

	#calculate the length of root trunk
	#remove the sample which have 0 muatation
	print(all_data)
	for item in dic.keys():
		if len(dic[item])==0:
			sample_shu=sample_shu-1
			print(item+' have '+str(len(dic[item]))+' mutation,so remove this sample')
		else:
			print(item+' have mutations: '+str(len(dic[item])))
	for key in list(dic.keys()):
		if not dic.get(key):
			del dic[key]
			all_data.remove(key)
	if len(all_data)>=2:
		all_driver_mutation={}
		for item in dic.keys():
			if len(search_gene(dic[item],'FALSE'))>1:
				for info in search_gene(dic[item],'TRUE')[:-1]:#remove count_map_driver_gene
					split=info.split(':')
					if split[0] not in all_driver_mutation.keys():
						all_driver_mutation.update({split[0]:[split[1]]})
					else:
						if split[1] not in all_driver_mutation[split[0]]:
							all_driver_mutation[split[0]].append(split[1])
		un_uniq_gene=0
		for key in all_driver_mutation.keys():
			if len(all_driver_mutation[key])>1:
				print(key+' have '+str(len(all_driver_mutation[key]))+' uniq mutations: '+str(all_driver_mutation[key]))
			else:
				un_uniq_gene=un_uniq_gene+1
		if un_uniq_gene==len(all_driver_mutation.keys()):
			print('every gene refer to a uniq gene')
		#........................................................................
		i=0
		for item in dic.keys():
			i=i+1
			if i==1:
				all_jiaoji=dic[item]
			all_jiaoji = [val for val in all_jiaoji if val in dic[item]]
		root_node=str(len(all_jiaoji))+'['+str(search_gene(all_jiaoji,'False')[-1])+'-'
		for item in search_gene(all_jiaoji,'False')[:-1]:
			root_node=root_node+item+';'
		root_node=root_node[:-1]+']'
		for item in dic.keys():
			for item1 in all_jiaoji:
				dic[item].remove(item1)

		#...............................................body......................
		dic_baocun={}
		for item in all_data:
			dic_baocun.update({item:[]})
		for item in dic_baocun.keys():
			for item_data in dic[item]:
				dic_baocun[item].append(item_data)
		global all_DATA
		all_result_1=[]
		global branch_path
		branch_path=[]
		global path_string
		path_string=[]
		global final_path_count
		final_path_count=-1
		global DIC_result
		DIC_result={}
		DIC_count=0
		duli_count=0
		duli_state='false'
		chongfu_bianli=0
		while final_path_count!=len(branch_path) or len(path_string)!=len(branch_path):
			chongfu_bianli=chongfu_bianli+1
			dic={}
			for item in all_data:
				dic.update({item:[]})
			for item in dic.keys():
				for item_data in dic_baocun[item]:
					dic[item].append(item_data)
			all_DATA=[]
			all_DATA.append(all_data)
			result_1=[]
			path=','
			dic_result={}
			for all_item in all_DATA:
				dic1={}
				for item in all_item:
					dic1.update({item:[]})
				for item in dic1.keys():
					for item_data in dic[item]:
						dic1[item].append(item_data)
				#check if there are any-two sample that do not share any mutation?
				while len(all_item)>=2:
					sample_2={}
					count2=0
					cishu2=int(Cmn(len(all_item),2))
					while len(sample_2.keys())<cishu2:
						sel2=random.sample(all_item,2)
						original_data=sel2[0]+','+sel2[1]
						fan_original_data=sel2[1]+','+sel2[0]
						if fan_original_data not in sample_2.keys():
							sample_2.update({original_data:[]})
					for item_2 in sample_2.keys():
						item_2_split=item_2.split(',')
						sample_2_jiaoji=[val for val in dic1[item_2_split[0]] if val in dic1[item_2_split[1]]]
						if len(sample_2_jiaoji)!=0:
							count2=count2+1
					if count2>0:
						duli_count=duli_count+1	
						result_all=xunhuan(all_item)
						path_state='false'
						for result_path in result_all:
							same_path_count=0
							identify1=result_path.split(':')[0]
							for result_path_1 in result_all:
								if identify1 == result_path_1.split(':')[0]:
									same_path_count=same_path_count+1
							if same_path_count>=2:
								path_state='true'
						branch_path_iter=[]
						for path_item in branch_path:
							if path_item not in branch_path_iter:
								branch_path_iter.append(path_item)
							else:
								print('branch_path error')
						#print('path_state: '+path_state)
						for result_all_1 in result_all:
							if path_state=='false':
								path_merge=path+';'+result_all_1.split(':')[0]
							else:
								path_merge=path+';'+result_all_1
							#print('path_merge: '+path_merge)
							if path_merge[0]==',' and chongfu_bianli==1:
								if path_merge[2:] not in branch_path:
									branch_path.append(path_merge[2:])
							else:
								if path_merge[0]==',':
									path_merge=path_merge[2:]
								path_all=len(branch_path_iter)
								path_count=0
								while path_count<path_all:
									if len(branch_path_iter[path_count].split(';'))==len(path_merge.split(';'))-1 and branch_path_iter[path_count] ==path_merge[:len(branch_path_iter[path_count])] :
										if  branch_path_iter[path_count] == path_merge[:len(branch_path_iter[path_count])]:
											if branch_path_iter[path_count] in branch_path:
												branch_path.remove(branch_path_iter[path_count])
										if path_merge not in branch_path:
											branch_path.append(path_merge)
									path_count=path_count+1			
		
						for result_all_1 in result_all:
							if path_state=='false':
								if path[0]==',':
									#print('ok1, in branch_path')
									path_merge=path+';'+result_all_1.split(':')[0]
									path_merge=path_merge[2:]
									path_copy_state='FALSE'
									for path_copy1 in branch_path:
										if path_copy1[-3:]!='end' and path_merge==path_copy1[:len(path_merge)] :
											path=path_merge
											path_copy_state='TRUE'
											break
									if path_copy_state=='TRUE':
										break
								else:
									path_merge=path+';'+result_all_1.split(':')[0]
									path_copy_state='FALSE'
									for path_copy1 in branch_path:
										if path_copy1[-3:]!='end' and path_merge==path_copy1[:len(path_merge)] :
											path=path_merge
											path_copy_state='TRUE'
											break
									if path_copy_state=='TRUE':
										break
							else:
								if path[0]==',':
									#print('ok3, in branch_path')
									path_merge=path+';'+result_all_1
									path_merge=path_merge[2:]
									path_copy_state='FALSE'
									for path_copy1 in branch_path:
										if path_copy1[-3:]!='end' and path_merge==path_copy1[:len(path_merge)] :
											path=path_merge
											path_copy_state='TRUE'
											break
									if path_copy_state=='TRUE':
										break
								else:
									path_merge=path+';'+result_all_1
									path_copy_state='FALSE'
									for path_copy1 in branch_path:
										if path_copy1[-3:]!='end' and path_merge==path_copy1[:len(path_merge)] :
											path=path_merge
											path_copy_state='TRUE'
											break
									if path_copy_state=='TRUE':
										break
						result=result_all_1
						result_split=result.split(':')
						result_split=result_split[1].split(',')
						first='.'
						#print('max_pair:')
						for item_split in result_split:
							first=first+','+item_split
						first=first[2:]
						return_2=[]
						for item_split in all_item:
							if item_split not in result_split:
								return_2.append(item_split)
						return_2.sort()
						#print(return_2)
						return_2_zuhe=[]
						while len(return_2_zuhe) <fact(len(return_2)):
							sel=random.sample(return_2,len(return_2))
							da=','
							for Item in sel:
								da=da+','+Item
							if da[2:] not in return_2_zuhe:
								return_2_zuhe.append(da[2:])
						#print(return_2_zuhe)
						second='.'
						for item3 in return_2_zuhe:
							if item3 in first_3.keys():
								second=item3
								break

						if int(first_3[first].split('[')[0])!=len(first_4[first]):
							print('false')
						result_item_1=first+':'+str(len(first_4[first]))+'['+str(search_gene(first_4[first],'False')[-1])+'-'
						for item in search_gene(first_4[first],'False')[:-1]:
							result_item_1=result_item_1+item+';'
						result_item_1=result_item_1[:-1]+']'
						result_1.append(result_item_1)
						dic_result.update({first:first_4[first]})
						result_item_2=second+':'+str(len(first_4[second]))+'['+str(search_gene(first_4[second],'False')[-1])+'-'
						for item in search_gene(first_4[second],'False')[:-1]:
							result_item_2=result_item_2+item+';'
						result_item_2=result_item_2[:-1]+']'
						result_1.append(result_item_2)
						dic_result.update({second:first_4[second]})
						if int(result[0:1])<2:
							result_split=result.split(':')
							result_split_split=result_split[1].split(',')
							all_item=[]
							for item in result_split_split:
								all_item.append(item)
							for key in list(dic1.keys()):
								if key not in result_split_split:
									del dic1[key]
							#filter key which value=null
							for key in list(dic1.keys()):
								if not dic1.get(key):
									del dic1[key]
									all_item.remove(key)
							#print(all_item)
							i=0
							for item_1 in dic1.keys():
								i=i+1
								if i==1:
									jiaoji=dic1[item_1]
								jiaoji=[val for val in jiaoji if val in dic1[item_1]]

							for item in dic1.keys():
								for item1 in jiaoji:
									dic1[item].remove(item1)

							for item in dic1.keys():
								for item1 in jiaoji:
									dic[item].remove(item1)
						else:
							group1=[]
							group2=[]
							result_split=result.split(':')
							result_split_split=result_split[1].split(',')
							for item in result_split_split:
								group1.append(item)
							group1.sort()
							i=0
							for item_1 in group1:
								i=i+1
								if i==1:
									jiaoji=dic[item_1]
								jiaoji=[val for val in jiaoji if val in dic[item_1]]

							for item in group1:
								for item1 in jiaoji:
									dic[item].remove(item1)

							for item in all_item:
								if item not in result_split_split:
									group2.append(item)
							group2.sort()
							i=0
							for item_1 in group2:
								i=i+1
								if i==1:
									jiaoji=dic[item_1]
								jiaoji=[val for val in jiaoji if val in dic[item_1]]

							for item in group2:
								for item1 in jiaoji:
									dic[item].remove(item1)
			
							all_DATA.append(group1)
							all_DATA.append(group2)

							break
					else:
						for item in all_item:
							result_item =item+':'+str(len(dic1[item]))+'['+str(search_gene(dic1[item],'False')[-1])+'-'
							for item0 in search_gene(dic1[item],'False')[:-1]:
								result_item=result_item+item0+';'
							result_item=result_item[:-1]+']'
							#print(result_item)
							result_1.append(result_item)
							dic_result.update({item:dic1[item]})
							if duli_count==0:
								result_item=result_item+';end'
								branch_path.append(result_item)
						if duli_count==0:
							all_result_1.append(result_1)
							duli_state='true'
						#print(sample_name)
						break
			if duli_count>0:
				path=path+';end'
				if path==',;end':
					print('error: '+path)
				for path_2 in branch_path:
					if path[:-4]==path_2:
						branch_path.remove(path_2)
				if path not in branch_path:
					branch_path.append(path)
			for i in branch_path:
				if branch_path.count(i)>1:
					print('branch_path have chongfu item: '+i)
			final_path_count=0
			for path_data in branch_path:
				if path_data.split(';')[-1]=='end' or path_data[0].isalpha():
					final_path_count=final_path_count+1

			if duli_state=='false':
				if path not in path_string and path[-3:]=='end':
					path_string.append(path)

					DIC_result.update({DIC_count:dic_result})
					DIC_count=DIC_count+1
					all_result_1.append(result_1)
			else:
				for path_1 in branch_path:
					if path_1 not in path_string and path_1[-3:]=='end':
						path_string.append(path_1)
					else:
						print('The tree structure that completely branches from the root node is wrong!')
				DIC_result.update({DIC_count:dic_result})
				DIC_count=DIC_count+1
		print('branch_path count: '+str(len(branch_path)))
		print('path_string: '+str(len(path_string)))
		print('all_result: '+str(len(all_result_1)))
		print('DIC_result: '+str(len(DIC_result.keys())))
		for item in all_result_1:
			item_copy = item[:]

			for data in item_copy:
				if int(data.split(':')[1].split('[')[0])==0 and len(data.split(':')[0].split(','))>1:
					item.remove(data)
		all_result_weight=[]
		for item in all_result_1:
			weight=0

			for data in item:
				length=len(data.split(':')[0].split(','))
				data=data.split(':')[1].split('[')[0]
				if length>=2:
					weight=weight+length*int(data)
			all_result_weight.append(weight)

		i=0
		max_weight_data=0
		while i <len(all_result_weight):
			if i==0:
				max_weight_data=all_result_weight[i]
			if all_result_weight[i]>max_weight_data:
				max_weight_data=all_result_weight[i]
			i=i+1
		result_max_weight=[]
		if max_weight_data==0:
			DIC_result_max_weight={}
			for item in all_result_1:
				result_max_weight.append(item)
			DIC_result_max_weight.update({0:DIC_result[0]})
		else:
			DIC_result_max_weight={}
			i=0
			i1=0
			while i <len(all_result_weight):
				if all_result_weight[i]==max_weight_data:
					all_result_1[i].sort()
					if all_result_1[i] not in result_max_weight:
						result_max_weight.append(all_result_1[i])
						DIC_result_max_weight.update({i1:DIC_result[i]})

						i1=i1+1
				i=i+1

		tree_count=0
		for result_1 in result_max_weight:
			dic_result=DIC_result_max_weight[tree_count]
			dic_final={}
			ll=[]#header
			ll.append('mutation_id')
			for item in dic_result.keys():
				if len(dic_result[item])!=0 or len(item.split(','))==1:
					ll.append(item)#all branch
			for item in ll[1:] :
				for mut in dic_result[item]:
					if mut not in dic_final.keys():
						LL=[]
						LL.append(mut)
						for item1 in ll[1:]:
							if mut in dic_result[item1]:
								LL.append(1)
							else:
								LL.append(0)
						driver_string='.'
						for driver_item in search_gene([mut],'FALSE')[:-1]:
							driver_string=driver_string+driver_item
						if len(driver_string)==1:
							LL.append('.')
						else:
							LL.append(driver_string[1:])
						dic_final.update({mut:LL})
			#print(dic_final)
			ll.append('driver_gene')
			result_2=[]
			number=1
			#print(str(sample_shu))
			while number<=sample_shu:
				for item in result_1:
					if number==1:
						if len(item.split(':')[0].split(','))==number:
							#print(item.split(':')[0].split(','))
							result_2.append(item)
					else:
						if len(item.split(':')[0].split(','))==number:						
							new_item=','
							if int(item.split(':')[1].split('[')[0])!=0:
								for item1 in item.split(':')[0].split(','):
									for item2 in result_2:
										num=0
										while num<len(item2):
											if item2[num]=='(':
												num=num+1
											else:
												break
										if item2[len(str(item1))+num]==':':
											if item1 ==item2[num:len(str(item1))+num]: 
												new_item=new_item+','+item2
												result_2.remove(item2)
								new_item='('+new_item[2:]+'):'
								for item3 in item.split(':')[1:]:
									new_item=new_item+item3+':'
								#print(new_item[:-1]) 
								result_2.append(new_item[:-1])
					#print(result_2)
				number=number+1
			#print(result_2)
			final_result=','
			for item in result_2:
				final_result=final_result+','+item
			final_result='('+final_result[2:]+'):'+root_node
			print('tree_structure:')
			print(final_result)
			if tree_count==0:
				first_final=final_result
			tree_structure=re.sub(u"\\[.*?]", "", final_result)
			tree_structure=tree_structure+';'
			#print(tree_structure)
			#.......................................plot.........................................
			if tree_count==0 or (tree_count>0 and final_result!=first_final):
				out=open(file_path_out+os.path.sep+patient+'_info_'+str(tree_count)+'.txt','w')
				for item in ll:
					out.write(str(item))
					out.write('\t')
				out.write('\n')
				for item in dic_final.keys():
					count=0
					for data in dic_final[item][1:-1]:
						count=count+int(data)
					if count>=2:
						out.write(item)
						out.write('\t')
						for data in dic_final[item][1:]:
							out.write(str(data))
							out.write('\t')
						out.write('\n')
				out.close()
				tree = Phylo.read(StringIO(tree_structure), "newick")
				#tree.ladderize()# Flip branches so deeper clades are displayed at top
				tree.rooted = True
				tree.name=patient
				tree = tree.as_phyloxml()
				tree.root.branch_labels=root_node
				#global tree_label
				if ('[' in final_result[1:].split(':')[-1]) and (']' in final_result[1:].split(':')[-1]):
					tree_label=final_result[1:].split(':')[:-1]
				else:
					tree_label=final_result[1:].split(':')[:-2]
				tree_label_1=[str(i) for i in tree_label]
				tree_labels = ':'.join(tree_label_1)
				tree_label_split=tree_labels[:-1]
				global new_tree_labels
				new_tree_labels=split_tree(tree_label_split)
				original_count=len(new_tree_labels)
				original_count1=original_count
				num=0
				for item in new_tree_labels:
					if num <=original_count-1:
						tree.clade[num].branch_labels=item.split(':')[-1]
						if ('(' in item) or (')' in item):
							new_item_1=[str(i) for i in item.split(':')[:-1]]
							new_item_split=':'.join(new_item_1)
							new_item_split=new_item_split[1:-1]
							new_items=split_tree(new_item_split)
							index=0
							for branch in new_items:
								branch=str(num)+'*'+str(index)+'*'+branch
								new_tree_labels.append(branch)
								index=index+1
						original_count1=original_count1-1
					else:
						branch=item.split('*')[-1]
						branch_count=item.split('*')[:-1]
						#print(branch_count)
						branch_count_merge=','
						for index in branch_count:
							branch_count_merge=branch_count_merge+'*'+index
						branch_count_merge=branch_count_merge[2:]
						str_count='tree.clade['
						for ITEM in branch_count:
							str_count=str_count+ITEM+','
						str_count=str_count[:-1]+'].branch_labels = '+'"'+item.split(':')[-1]+'"'
						#print(str_count)
						exec(str_count)
						#print(tree)
						if ('(' in item) or (')' in item):
							new_item_1=[str(i) for i in item.split(':')[:-1]]
							new_item_split=':'.join(new_item_1)
							new_item_split=new_item_split.split('*')[-1]
							new_item_split=new_item_split[1:-1]
							new_items=split_tree(new_item_split)
							index=0
							for branch1 in new_items:
								branch1=str(branch_count_merge)+'*'+str(index)+'*'+branch1
								new_tree_labels.append(branch1)
								index=index+1
					num=num+1
				matplotlib.rc('font', size=6)
				plt.rcParams['lines.linewidth'] = 0.7
				Phylo.draw(tree,do_show=False,branch_labels=lambda c: c.branch_labels)	
				#Phylo.draw(tree,do_show=False,branch_labels=lambda c: int(c.branch_length))
				plt.savefig(file_path_out+os.path.sep+patient+'_tree_'+str(tree_count)+'.pdf',dpi=500)
				plt.close()
				tree_count=tree_count+1
		print(patient+' ok!')
	else:
		print(patient+' have less patients!')
