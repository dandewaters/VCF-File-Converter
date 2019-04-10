# VCF-File-Converter
Converts files in Variant-Calling Format (VCF) to Heterozygous Allele Depth Format (HAD), the format required for GBS2Ploidy. Also converts VCF files to the format for COLONY software.


For converting to colony format:
 
I have allowed for this program to separate data for offspring and (potential) maternal specimens. If you want the program to automatically make separate maternal and offspring files, you must use our naming convention. We label our maternal individuals with their group name followed by their number in the group. For example, these sample names ("abc1", "abc2" , "abc3") are for maternal specimens from project "abc", with numbers 1, 2, and 3. The offspring will have the maternal individual's name followed by the letter "e" (for egg) and then followed by that individual's number in the group. For example, if the maternal individual's name is "abc1", their offspring will be named "abc1e1", "abc1e2", "abc1e3", etc.
 
 What's important is that the offspring individuals' names must end in a number, followed by the letter "e", followed by another number. I used the regular expressions library to search through the individuals names in the vcf files with the string "\de\d*?". If that search returns a match, the program writes that individual's information into the offspring file. otherwise, it goes into the maternal file.
