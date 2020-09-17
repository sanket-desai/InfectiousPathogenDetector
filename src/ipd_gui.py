#!/usr/bin/python3
'''

Author	  : Sonal Rashmi
Date		: 16/07/2020
Description : IPD GUI enabled with two layer frame shift (Multisample/Singlesample and Single-end/Paired-end)

'''

from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import os.path
from commandlineparser import *
from ipdshortread import *
from ipdlongread import *
import subprocess

#All methods describing widgets and their functions
def outputdirectory_browse_button():
	folder_name = filedialog.askdirectory(initialdir = "/")
	ouputdirectory_temp.set(folder_name)
		
def inputdirectory_browse_button():
	folder_name = filedialog.askdirectory(initialdir = "/")
	inputdirectory_temp.set(folder_name)
	
def label(master,text,grid_row,grid_col):
	lbl = Label(master=master,text=text)
	lbl.grid(row=grid_row, column=grid_col)

def text_box(master,textvariable,grid_row,grid_col):
	txt = Entry(master=master,textvariable=textvariable)
	txt.grid(row=grid_row, column=grid_col)

def button(master,txt,cmd,grid_row,grid_col):
	button = Button(master=master,text=txt, command=cmd)
	button.grid(row=grid_row, column=grid_col, sticky = W+E)
	
def radioButton(master,radio_button_label_value_dict,grid_row,grid_col,rad_variable):
	for (text, value) in radio_button_label_value_dict.items(): 
		rad = Radiobutton(master=master,text = text, variable = rad_variable, value = value)
		rad.grid(row=grid_row, column=grid_col)
		grid_col=int(grid_col)+1
	rad_variable.set(value)
	
def scroll_bar(master,int_from,int_to,grid_row,grid_col,spin_variable,default=None,width=None):
	if default:
		spin_variable.set(default)
	else:
		spin_variable.set(int_from)
	if width: 
		spin = Spinbox(master=master, from_=int_from, to=int_to, width=width,textvariable=spin_variable)
	else:
		spin = Spinbox(master=master, from_=int_from, to=int_to, textvariable=spin_variable)
	spin.grid(row=grid_row, column=grid_col)	

def input_fastq_file_browse_button():
	file = filedialog.askopenfilenames(initialdir = "/", filetypes = [('Fastq Files', '*.fq'),('Fastq Files', '*.fq.gz'),('Fastq Files', '*.fastq'),('Fastq Files', '*.fastq.gz')]) 
	fastq_file_input_temp.set(file)
	
def input_fastq_file_browse_button_2():
	file = filedialog.askopenfilenames(initialdir = "/", filetypes = [('Fastq Files', '*.fq'),('Fastq Files', '*.fq.gz'),('Fastq Files', '*.fastq'),('Fastq Files', '*.fastq.gz')]) 
	fastq_file_input_temp_2.set(file)	

def input_all_file_browse_button():
	file = filedialog.askopenfilenames(initialdir = "/", filetypes = [("all files","*.*")]) 
	file_input_temp.set(file)

#Description of all the 4 frames and their functions

def singlessframe():
	global fastq_file_input_temp
	
	fastq_file_input_temp = StringVar()
	
	single_ps_frame.grid_forget()
	single_ss_frame.grid_propagate(1)
	single_ss_frame.grid(row=7, column=1)	
	
	lbl8_2_1 = label(master=single_ss_frame,text="Fastq File",grid_row=7, grid_col=0)
	button8_2_1 = button(master=single_ss_frame, txt="Browse", cmd=input_fastq_file_browse_button,grid_row=7, grid_col=2)
	txt8_2_1 = text_box(master=single_ss_frame,textvariable=fastq_file_input_temp,grid_row=7, grid_col=1)

	return(fastq_file_input_temp)

def singlepsframe():
	global fastq_file_input_temp
	global fastq_file_input_temp_2
	
	fastq_file_input_temp = StringVar()
	fastq_file_input_temp_2 = StringVar()
	
	single_ss_frame.grid_forget()
	single_ps_frame.grid_propagate(1)
	single_ps_frame.grid(row=7, column=1)
	
	lbl8_1_1 = label(master=single_ps_frame,text="Fastq File",grid_row=7, grid_col=1)
	button8_1_1 = button(master=single_ps_frame, txt="Browse", cmd=input_fastq_file_browse_button,grid_row=7, grid_col=3)
	txt8_1_1 = text_box(master=single_ps_frame,textvariable=fastq_file_input_temp,grid_row=7, grid_col=2)

	lbl8_1_2 = label(master=single_ps_frame,text="Fastq File",grid_row=8, grid_col=1)
	button8_1_2 = button(master=single_ps_frame, txt="Browse", cmd=input_fastq_file_browse_button_2,grid_row=8, grid_col=3)
	txt8_1_2 = text_box(master=single_ps_frame,textvariable=fastq_file_input_temp_2,grid_row=8, grid_col=2)
	
	return ([fastq_file_input_temp,fastq_file_input_temp_2])
	
def multiframe():
	global file_input_temp
	
	file_input_temp = StringVar()

	single_frame.grid_forget()
	multi_frame.grid_propagate(1)
	multi_frame.grid(row=6, column=4)
	
	lbl7_1_1 = label(master=multi_frame,text="Sample File",grid_row=6, grid_col=1)
	button7_1_1 = button(master=multi_frame, txt="Browse", cmd=input_all_file_browse_button,grid_row=6, grid_col=3)
	txt7_1_1 = text_box(master=multi_frame,textvariable=file_input_temp,grid_row=6, grid_col=2)



def singleframe():
	global data_type
	data_type=StringVar()
	
	multi_frame.grid_forget()
	single_frame.grid_propagate(1)
	single_frame.grid(row=6, column=6)
	
	lbl7_2_1 = label(master=single_frame,text="Data Type",grid_row=6, grid_col=0)
	Radiobutton(master=single_frame, text="Single-ended Samples", variable=data_type, value='se', command=singlessframe).grid(column=1, row=6, sticky = W+E)
	Radiobutton(master=single_frame, text="Paired-ended Sample", variable=data_type, value='pe', command=singlepsframe).grid(column=2, row=6, sticky = W+E)

#Descriptions of all the error checks and warnings  

def show_warning_message_box(warning_message):
	messagebox.showwarning("Warning",warning_message)

def check_file_exists(Variable):
	if os.path.isfile(Variable):
		return Variable
	else:
		show_warning_message_box("File " +str(Variable) + " doesn't exist")
		
def check_directory_exists(Variable):
	if os.path.exists(Variable):
		return Variable
	else:
		show_warning_message_box("Directory " +str(Variable) + " doesn't exist")

def check_variable_exists(Variable):
	if Variable:
		return Variable
	else:
		show_warning_message_box(str(Variable) + " Empty Variable")


def reset_all():
	ouputdirectory_temp.set("")
	project_name_temp.set("")
	threads_temp = 4
	run_mode_temp.set(None)
	sequencing_type_temp.set(None)
	data_type.set(None)
	fastq_file_input_temp.set("")
	fastq_file_input_temp_2.set("")
	file_input_temp.set("")
	
#Input checks and creation of imaps objects 
		
def input_checks(master):
	OuputDirectory = None
	ProjectName = None
	Threads = None
	RunMode = None
	DataType	= None
	FastqFileInput = None
	SampleInput = None
	FastqFileInput2 = None
	SequencingType = None
	MolecularType = None
	CheckReport = None
	RunMode=check_variable_exists(run_mode_temp.get())
	OuputDirectory=check_directory_exists(ouputdirectory_temp.get())
	ProjectName=check_variable_exists(project_name_temp.get())
	SequencingType=check_variable_exists(sequencing_type_temp.get())
	MolecularType=check_variable_exists(molecular_type_temp.get())
	Threads=check_variable_exists(threads_temp.get())
	CheckReport=check_variable_exists(checkreport_temp.get())
	if RunMode == 'l':
		SampleInput=check_file_exists(file_input_temp.get().split(',')[0].split('(')[1].replace("'", ""))
		return_list=["Multisample", OuputDirectory, ProjectName, MolecularType, Threads, SampleInput, FastqFileInput, FastqFileInput2]
		Parser=CommandLineInterfaceParser(return_list)
		maplist=Parser.map_multiple_sample()
	elif RunMode == 's':
		DataType=data_type.get()
		if DataType == "se":
			FastqFileInput=check_file_exists(fastq_file_input_temp.get().split(',')[0].split('(')[1].replace("'", ""))
			return_list=["SingleEndSample", OuputDirectory, ProjectName, MolecularType, Threads, SampleInput, FastqFileInput, FastqFileInput2]
			Parser=CommandLineInterfaceParser(return_list)
			maplist=Parser.map_single_end_sample()
		elif DataType == "pe":
			if SequencingType == 'long':
				show_warning_message_box("Option not Available")
			else:	
				FastqFileInput=check_file_exists(fastq_file_input_temp.get().split(',')[0].split('(')[1].replace("'", ""))
				FastqFileInput2=check_file_exists(fastq_file_input_temp_2.get().split(',')[0].split('(')[1].replace("'", ""))
				return_list=["PairedEndSample", OuputDirectory, ProjectName, MolecularType, Threads, SampleInput, FastqFileInput, FastqFileInput2]
				Parser=CommandLineInterfaceParser(return_list)
				maplist=Parser.map_paired_sample()			
		else:
			show_warning_message_box("Invalid Option for DataType")
			sys.exit(0)
	else:
		show_warning_message_box("Invalid Option for RunMode")
		sys.exit(0)
	boo=0	
	if len(maplist) > 0:
		if SequencingType == 'long':
			for inmap in maplist: 					
				try:
					i=IPDLongRead(inmap)
					if CheckReport == "1":
						#print(CheckReport)
						cmd="python3 cov2reportgenerator.py -dir "+OuputDirectory
						cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
						cprocess.check_returncode()
					boo=1
				except TypeError:
					show_warning_message_box("Error in code !!!!")
		elif SequencingType == 'short':
			for inmap in maplist: 					
				try:
					i=IPDShortRead(inmap)
					if CheckReport == "1":
						#print(CheckReport)
						cmd="python3 cov2reportgenerator.py -dir "+OuputDirectory
						cprocess=subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)
						cprocess.check_returncode()
					boo=1
				except TypeError:
					show_warning_message_box("Error in code !!!!")
	if boo:
		a=None
		if messagebox.askokcancel("Success","Successfully Completed. Do you want to reset?",default="cancel"):
			reset_all()
	#print(return_list)	

def main():	
	
	root = Tk()
	root.title("IPD")
	root.geometry("1000x300")
	root.resizable(width=FALSE, height=FALSE)


	global ouputdirectory_temp
	global data_type_temp
	global project_name_temp
	global threads_temp
	global run_mode_temp
	global data_type
	global sequencing_type_temp
	global molecular_type_temp
	global multi_frame 
	global single_frame 
	global single_ss_frame 
	global single_ps_frame
	global fastq_file_input_temp
	global fastq_file_input_temp_2
	global file_input_temp
	global checkreport_temp

	ouputdirectory_temp = StringVar()
	project_name_temp=StringVar()
	threads_temp=IntVar()
	sequencing_type_temp=StringVar()
	molecular_type_temp=StringVar()
	run_mode_temp=StringVar()
	fastq_file_input_temp=StringVar()
	fastq_file_input_temp_2=StringVar()
	file_input_temp=StringVar()
	data_type=StringVar()
	checkreport_temp=IntVar()
	
	menu=Menu(root) #Creates a menu button
	root.config(menu=menu)
	#Creating the file menu
	file_item=Menu(menu)
	menu.add_cascade(label='File', menu=file_item)
	menu.add_command(label='About')
	menu.add_command(label='Exit', command=root.quit)


	sequencing_type_dict={"Short Read":'short', "Long Read":'long'}
	molecular_type_dict={"DNA":'DNA', "RNA":'RNA'}

	lbl1 = label(master=root,text="Output Directory",grid_row=0, grid_col=0)
	button1 = button(master=root,txt="Browse", cmd=outputdirectory_browse_button,grid_row=0, grid_col=2)
	txt1 = text_box(master=root,textvariable=ouputdirectory_temp,grid_row=0, grid_col=1)

	lbl2 = label(master=root,text="Project Name",grid_row=1, grid_col=0)
	txt2 = text_box(master=root,textvariable=project_name_temp,grid_row=1, grid_col=1)

	lbl3 = label(master=root,text="Threads",grid_row=2, grid_col=0)
	scroll_bar1=scroll_bar(master=root,int_from=1,int_to=100,grid_row=2,grid_col=1,spin_variable=threads_temp,default=4)

	lbl4 = label(master=root,text="Sequencing Type",grid_row=3, grid_col=0)
	rad_button1=radioButton(master=root,radio_button_label_value_dict=sequencing_type_dict,grid_row=3, grid_col=1,rad_variable=sequencing_type_temp)

	lbl5 = label(master=root,text="Molecular Type",grid_row=4, grid_col=0)
	rad_button2=radioButton(master=root,radio_button_label_value_dict=molecular_type_dict,grid_row=4, grid_col=1,rad_variable=molecular_type_temp)

	lbl6 = label(master=root,text="Run Mode",grid_row=5, grid_col=0)
	Radiobutton(master=root, text="Multi Samples", variable=run_mode_temp, value='l', command=multiframe).grid(column=1, row=5, sticky = W+E)
	Radiobutton(master=root, text="Single Sample", variable=run_mode_temp, value='s', command=singleframe).grid(column=2, row=5, sticky = W+E)

	multi_frame = Frame(master=root,  width=1, height=2, bd=1, relief=SUNKEN)
	single_frame = Frame(master=root,  width=1,height=2, bd=1, relief=SUNKEN)

	single_ss_frame = Frame(master=single_frame, width=1, height=2, bd=1, relief=SUNKEN)
	single_ps_frame = Frame(master=single_frame, width=1, height=2, bd=1, relief=SUNKEN)
	
	lbl9 = label(master=root,text="Generate SARS-CoV2 Report",grid_row=9, grid_col=0)
	checkbox1 = Checkbutton(master=root, variable=checkreport_temp, onvalue="1", offvalue="0")
	checkbox1.grid(row=9, column=1, sticky = W+E)
	
	button2=button(master=root,txt="Submit",cmd=lambda: input_checks(master=root),grid_row=10, grid_col=1) 

	mainloop()
	
if __name__ =="__main__":
	main()
	
