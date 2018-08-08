import os

# Specimen class, used for the purpose of organizing the data in these files
# Consists of a string to repesent the name of the specimen (E.g. 10a), and
# a list of strings containing the loci information.
class Specimen():
	
	# Initializes specimen object
	def __init__(self, name):
		self.m_name = name
		self.m_loci_list = []

	# SETTERS
	# Sets name
	def set_name(self, name):
		self.m_name = name
	# Sets loci list
	def set_loci_list(self, loci_list):
		self.m_loci_list = loci_list
	# Sets loci at a specific index in the loci list
	def set_loci_at(self, loci, index):
		self.m_loci_list[index] = loci

	# GETTERS
	# Returns specimen name
	def get_name(self):
		return(self.m_name)
	# Returns specimen loci list
	def get_loci_list(self):
		return(self.m_loci_list)
	# Returns specimen loci at a specific index
	def get_loci_at(self, index):
		return(self.m_loci_list[index])	

	# OTHER
	# Adds loci to the loci list
	def add_loci(self, loci):
		self.m_loci_list.append(loci)

	# (For debugging purposes)
	def print_data(self):
		print("Name:", self.get_name())
		print("Loci:")
		for loci in self.get_loci_list():
			print("\t", loci)


# Opens VCF files in the VCF folder, makes a list of the lines in the file
# and removes the commands at the top of the file and the 9 leftmost columns
def openVCF(VCFFileName):

	print("\nNow reading:", VCFFileName)
	VCFFile = open("C:\\Users\\BurnsLab\\Desktop\\VCF Converter Tool\\VCF Files\\" + VCFFileName, "r")
	# Creates a list of strings of file contents, line by line
	#(Clears list when reading a new file)
	VCFFileContents = []
	VCFFileContents = VCFFile.readlines()
	
	# Removes commands/file info at the top
	while(VCFFileContents[0].startswith("##")):
		VCFFileContents.remove(VCFFileContents[0])
	
	# Splits line into columns and removes the first 9 columns in each line
	for line in range(0, len(VCFFileContents)):
		# Strips off "\n" at the end of each line
		VCFFileContents[line] = VCFFileContents[line].strip()
		# Splits list into 2-dimensional list
		VCFFileContents[line] = VCFFileContents[line].split("\t")
		# Slices off commands at top of file
		VCFFileContents[line] = VCFFileContents[line][9:]

	VCFFile.close()
	return(VCFFileContents)


# Replaces homozygous loci with ["NA", "NA"] and heterozygous with nonzero numbers
def analyzeLociDataList(intermediateList):
	
	print("Searching for heterozygotes...")
	# Searches for elements of the list (skipping the specimen names)
	# that do not start with "./." or "0/0"
	for row in range(1, len(intermediateList)):
		for column in range(0, len(intermediateList[row])):
			if(intermediateList[row][column][0] == intermediateList[row][column][2]):
				intermediateList[row][column] = ["NA","NA"]
			else:
				# Splits loci data by ":"
				intermediateList[row][column] = intermediateList[row][column].split(":")
				# Removes first 2 elements of the list made from splitting
				intermediateList[row][column] = intermediateList[row][column][2]
				# Splits remaining numbers, adds nonzero numbers to a temp list
				intermediateList[row][column] = intermediateList[row][column].split(",")
				nonZeroNums = []
				for num in intermediateList[row][column]:
					if(int(num) > 0):
						nonZeroNums.append(num)
				
				# Keeps appropriate heterozygote information in list
				if(len(nonZeroNums) < 2):
					# Just in case pyrad fasely labels a heterozygote
					intermediateList[row][column] = ["NA","NA"]
				else:
					# Sets old loci data equal to the temp list
					intermediateList[row][column] = nonZeroNums
				
	return(intermediateList)


# Restructures the information from the list returned from analyzeLociDataList
# by making a list of specimen objects
def createSpecimen(lociList):

	print("Reorganizing heterozygote data...")
	specimenList = []
	
	# Iterates through first row of the loci list (row with specimen names)
	for column in range(0, len(lociList[0])):
		# Creates a specimen object from name in row 0 in lociList
		specimenList.append(Specimen(lociList[0][column]))
	# Removes name row in lociList
	lociList.remove(lociList[0])
	
	# Adds each loci to respective Specimen object
	for lociIndex in range(0, len(lociList)):
		for specIndex in range(0, len(specimenList)):
			specimenList[specIndex].add_loci(lociList[lociIndex][specIndex])

	# Returns list of specimen
	return(specimenList)


# Creates a .txt file, writes information from the specimen list to it
# in HetAlleleDepth format
def createTXTFile(fileName, specList):
	
	print("Writing data to HetAlleleDepth file")
	# Creates a respective name from the corresponding .VCF file
	fileName = fileName[:-4]
	fileName = "HetAlleleDepth" + fileName + ".txt"
	# Opens .txt file with respective name
	TXTFile = open("C:\\Users\\BurnsLab\\Desktop\\VCF Converter Tool\\Converted Files\\" + fileName, "w")
	
	# Writes names at the top of txt file
	line = ""
	for spec in specList:
		line += spec.get_name() + "\t" + spec.get_name() + "\t"
	line = line.strip()
	TXTFile.write(line + "\n")

	# Writes rest of data
	for lociIndex in range(len(specList[0].get_loci_list())):
		line = ""
		for spec in range(0, len(specList)):
			line += specList[spec].get_loci_at(lociIndex)[0] + "\t"
			line += specList[spec].get_loci_at(lociIndex)[1] + "\t"
		line = line.strip()
		TXTFile.write(line + "\n")
	
	TXTFile.close()



def main():
	
	print("Welcome to Dan's VCF file converting tool!")
	print("This program converts VCF files to Heterozygous Allele Depth Files.")
	print("Before we start converting, make sure every VCF file you want to convert is in the \"VCF Files\" folder")
	input("Hit Enter to start converting: ")
	
	#Creates list of files in the folder for VCF files to be converted
	VCFFileList = os.listdir("C:\\Users\\BurnsLab\\Desktop\\VCF Converter Tool\\VCF Files")
	print("Would you like to convert these files?\n")

	for VCFFileName in VCFFileList:
		print(VCFFileName)
	
	# Opens a VCF file in VCF folder, Removes unnecessary information,
	# restructures the information, and writes information to .txt file.
	# Repeats for each file in the folder
	readyToConvert = input("\nType \"y\" to continue, \"n\" to quit: ")
	if(readyToConvert == "y"):
		for VCFFileName in VCFFileList:
			lociDataList = openVCF(VCFFileName)
			analyzedList = analyzeLociDataList(lociDataList)
			specimenList = createSpecimen(analyzedList)
			createTXTFile(VCFFileName, specimenList)

	print("\n\nAll files have been converted successfully!")
	print("Navigate to: \"C:\\Users\\BurnsLab\\Desktop\\VCF Converter Tool\\Converted Files\" to access your converted files.")

main()
input("\n\nPress enter to quit")
