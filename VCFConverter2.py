import os

# Directories for input and output
VCFFOLDERDIR = "C:\\Users\\Dan\\OneDrive\\Dan School\\UMBC\\Lab\\File converter 2.0\\VCF Files\\"
CONVERTEDFILEFOLDER = "C:\\Users\\Dan\\OneDrive\\Dan School\\UMBC\\Lab\\File converter 2.0\\Converted Files\\"

# Specimen class, used for the purpose of organizing the data in these files
# Consists of a string to repesent the name of the specimen (E.g. 10a), and
# a list of strings containing the loci information.
class Specimen():
	
	# Initializes specimen object
	def __init__(self, name):
		self.m_name = name
		self.m_loci_list = []
		self.m_twoColumns = False
		self.m_allMissing = True

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

	# Sets boolean for if this specimen needs two columns in outputted file
	def set_two_columns(self):
		self.m_twoColumns = True
	# Sets boolean to signify that this specimen should not be included in outputted
	# file due to all loci being missing
	def set_all_missing(self):
		self.m_allMissing = False

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
	# Returns signal for necessity of two columns for this specimen
	def get_two_columns(self):
		return self.m_twoColumns
	# Returns signal for inclusion of this specimen in outputted file
	def get_all_missing(self):
		return self.m_allMissing

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
	VCFFile = open(VCFFOLDERDIR + VCFFileName, "r")
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



# Creates new TXT file, writes VCF data in HAD Format
def convertToHAD(FileContents, fileName, removeDoubleHets):
	
	print("Opening TXT file")
	# Creates a respective name from the corresponding .VCF file
	fileName = fileName[:-4]
	fileName = "HAD" + fileName + ".txt"
	# Opens .txt file with respective name
	TXTFile = open(CONVERTEDFILEFOLDER + fileName, "w")
	print("Writing to TXT file in HAD format...")
	
	# Writes names to file
	print("Writing names to output file...")
	for name in FileContents[0]:
		TXTFile.write(name + "\t" + name + "\t")
	TXTFile.write("\n")

	# Analyzes loci data and writes to txt file in HAD format
	for row in range(1, len(FileContents)):
		for column in range(0, len(FileContents[row])):
			if(FileContents[row][column][0] == FileContents[row][column][2]):
				TXTFile.write("NA" + "\t" + "NA" + "\t")
			else:
				# Splits loci data by ":"
				FileContents[row][column] = FileContents[row][column].split(":")
				# Removes first 2 elements of the list made from splitting
				FileContents[row][column] = FileContents[row][column][2]
				# Splits remaining numbers, adds nonzero numbers to a temp list
				FileContents[row][column] = FileContents[row][column].split(",")
				nonZeroNums = []
				for num in FileContents[row][column]:
					if(int(num) > 0):
						nonZeroNums.append(num)
				
				# Keeps appropriate heterozygote information in list
				if(len(nonZeroNums) < 2):
					# Just in case pyrad fasely labels a heterozygote
					TXTFile.write("NA" + "\t" + "NA" + "\t")
				elif (len(nonZeroNums) == 2):
					# Sets old loci data equal to the temp list
					TXTFile.write(nonZeroNums[0] + "\t" + nonZeroNums[1] + "\t")
				else:
					if(removeDoubleHets == True):
						TXTFile.write("NA" + "\t" + "NA" + "\t")
					else:
						nonZeroNums.sort()
						TXTFile.write(nonZeroNums[0] + "\t" + nonZeroNums[1] + "\t")

		# Starts new line
		TXTFile.write("\n")

	TXTFile.close()



# Rewrites data in appropriate format for colony
def analyzeLociDataList(intermediateList):
	print("Analyzing loci...")
	for row in range(1, len(intermediateList)):
		for column in range(0, len(intermediateList[row])):
			# Shorthand for first 3 characters of loci information
			read = intermediateList[row][column][0:3]
			# Case for heterozygotes
			if(read[0] != read[2]):
				intermediateList[row][column] = [int(read[0])+1, int(read[2])+1]
			# Homozygote case
			else:
				# Loci did not exceed threshold
				if(read == "./."):
					intermediateList[row][column] = [-9]
				else:
					intermediateList[row][column] = [int(read[0])+1]

	return(intermediateList)



# Restructures the information from the list returned from analyzeLociDataList
# by making a list of specimen objects
def createSpecimen(lociList):

	print("Reorganizing loci data...")
	specimenList = []
	
	# Iterates through first row of the loci list (row with specimen names)
	for column in range(0, len(lociList[0])):
		# Creates a specimen object from name in row 0 in lociList
		specimenList.append(Specimen(lociList[0][column]))
	# Removes name row in lociList
	lociList.remove(lociList[0])
	
	# Adds each loci to respective Specimen object
	for specIndex in range(0, len(specimenList)):
		for lociIndex in range(0, len(lociList)):
			if(len(lociList[lociIndex][specIndex]) > 1):
				specimenList[specIndex].set_two_columns()
				specimenList[specIndex].set_all_missing()
			elif(int(lociList[lociIndex][specIndex][0]) != -9):
				specimenList[specIndex].set_all_missing()
			specimenList[specIndex].add_loci(lociList[lociIndex][specIndex])

	# Returns list of specimen
	return(specimenList)



# Creates a .txt file, writes information from the specimen list to it
# in HetAlleleDepth format
def createTXTFile(fileName, specList):
	
	print("Writing data to colony file")
	# Creates a respective name from the corresponding .VCF file
	fileName = fileName[:-4]
	fileName = "COL" + fileName + ".txt"
	# Opens .txt file with respective name
	TXTFile = open(CONVERTEDFILEFOLDER + fileName, "w")
	
	# Writes names at the top of txt file
	line = ""
	for spec in specList:
		# Checks if specimen should be included in output file
		if(spec.get_all_missing() == False):
			line += spec.get_name() + "\t"
			# Writes name twice if two columns are necessary
			if(spec.get_two_columns()):
				line += spec.get_name() + "\t"
	line = line.strip()
	TXTFile.write(line + "\n")

	# Writes rest of data
	for lociIndex in range(len(specList[0].get_loci_list())):
		line = ""
		for spec in range(0, len(specList)):
			if(specList[spec].get_all_missing() == False):
				if(specList[spec].get_two_columns):
					if(len(specList[spec].get_loci_at(lociIndex)) > 1):
						line += str(specList[spec].get_loci_at(lociIndex)[0]) + "\t"
						line += str(specList[spec].get_loci_at(lociIndex)[1]) + "\t"
					else:
						line += (str(specList[spec].get_loci_at(lociIndex)[0]) + "\t") * 2
		line = line.strip()
		TXTFile.write(line + "\n")
	
	TXTFile.close()



def main():

	print("Welcome to Dan's VCF file converting tool 2.0!")
	print("This program converts VCF files to HAD or COLONY format.")
	print("Before we start converting, make sure every VCF file you want to convert is in the \"VCF Files\" folder")
	input("Hit Enter to start converting: ")

	#Creates list of files in the folder for VCF files to be converted
	VCFFileList = os.listdir(VCFFOLDERDIR)
	print("\n\nWould you like to convert these files?\n")

	for VCFFileName in VCFFileList:
		print(VCFFileName)
	
	# Opens a VCF file in VCF folder, Removes unnecessary information,
	# restructures the information, and writes information to .txt file.
	# Repeats for each file in the folder
	readyToConvert = input("\nType \"1\" to continue, anything else to quit: ")
	if(readyToConvert == "1"):
		print("\n\nWhat format do you want to convert these files to?")
		print("\t 1: Heterozygous Allele Depth Format")
		print("\t 2: COLONY Format")
		convertingFormat = input("\nType \"1\" or \"2\": ")
		
		# Converting to HAD format
		if(convertingFormat == "1"):
			print("\nDo you want to remove double heterozygotes?")
			print("\t 1: Yes")
			print("\t 2: No")
			remDubHet = input("\nType \"1\" or \"2\": ")

			# Begins converting VCF files in input folder to HAD format
			for VCFFileName in VCFFileList:
				VCFFileContents = openVCF(VCFFileName)

				if(remDubHet == "1"):
					convertToHAD(VCFFileContents, VCFFileName, True)
				else:
					convertToHAD(VCFFileContents, VCFFileName, False)
		
		# Converting to COLONY format
		elif(convertingFormat == "2"):
			for VCFFileName in VCFFileList:
				lociDataList = openVCF(VCFFileName)
				analyzedList = analyzeLociDataList(lociDataList)
				specimenList = createSpecimen(analyzedList)
				createTXTFile(VCFFileName, specimenList)

		print("\n\nAll files have been converted successfully!")
		print("Navigate to:", CONVERTEDFILEFOLDER, "to access your converted files.")

	print("Thank you for using my file converter tool!")

main()
