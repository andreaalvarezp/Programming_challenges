## ----------------------
## ASSIGNMENT 4
## Andrea Álvarez Pérez
## ----------------------

require 'bio'
require 'stringio'
input = ARGV

## BONUS:
## In order to find out whether the reciprocal best hits are actually between related genes, it would be nice to check
## if the gene functions (looking at GO Terms and KEGG) are actually similar or path-related.
## Also, in this kind of studies, we are interested in finding out the phylogenetic relationships of the sequences,
## maybe construction a phylogenetic tree of part of the sequences (the ones which give a reciprocal best hit) and some outliers.
## If tha analysis is correct, the outliers should be far from the match sequences in the tree.

# ref: "Choosing BLAST options for better detection of orthologs as reciprocal best hits." (Moreno-Hagelsieb G, Latimer K.)
# https://www.ncbi.nlm.nih.gov/pubmed/18042555
$evalue = 10**-6

puts "\n---------------------------------------------"
puts "ASSIGNMENT 4"
puts "Searching for orthologues"
puts "Andrea Álvarez Pérez"
puts "---------------------------------------------\n"

if input.length != 4
    puts "Incorrect number of arguments"
    puts "USAGE: $bioruby <filename1> <type filename1> <filename2> <type filename2>"
    puts "'Type' could be either 'nucl' or 'prot'"
    puts "---------------------------------------------\n";puts
    exit 1
end

file1 = input[0]
type_file1 = input[1]
file2 = input[2]
type_file2 = input[3]

if File.exists?("orthologues_output.txt")
    File.delete("orthologues_output.txt")
end
output = File.open("orthologues_output.txt", "a")
output.puts "----------------------------------------------------------"
output.puts "ASSIGMENT 4 - Searching for orthologues"
output.puts "Andrea Álvarez Pérez"
output.puts "----------------------------------------------------------";output.puts
output.puts "Orthologue pairs found between Arabidopsis and S. pombe";output.puts
output.puts "File 1 #{file1} ID\t----------\tFile 2 #{file2} ID"
output.puts "----------------------------------------------------------"

puts "\nFile 1: #{file1}";puts"\tType: #{type_file1}"
puts "File 2: #{file2}";puts"\tType: #{type_file2}"
puts "---------------------------------------------\n"

file1_db = file1.to_s.split(".")[0]+"_db"
file2_db = file2.to_s.split(".")[0]+"_db"

# We create the databases
system("mkdir Databases")
system("makeblastdb -in '#{file1}' -dbtype #{type_file1} -out ./Databases/#{file1_db.to_s}") 
system("makeblastdb -in '#{file2}' -dbtype #{type_file2} -out ./Databases/#{file2_db.to_s}")

# Depending on the seq type of the filesm we perform a type of blast
if type_file1 == 'nucl' and type_file2 == 'nucl'
    factory_f1 = Bio::Blast.local('blastn', "./Databases/#{file1_db.to_s}")
    factory_f2 = Bio::Blast.local('blastn', "./Databases/#{file2_db.to_s}")

elsif type_file1 == 'nucl' and type_file2 == 'prot' 
    factory_f1 = Bio::Blast.local('tblastn', "./Databases/#{file1_db.to_s}")
    factory_f2 = Bio::Blast.local('blastx', "./Databases/#{file2_db.to_s}")

elsif type_file1 == 'prot' and type_file2 == 'nucl'
    factory_f1 = Bio::Blast.local('blastx', "./Databases/#{file1_db.to_s}")
    factory_f2 = Bio::Blast.local('tblastn', "./Databases/#{file2_db.to_s}")

elsif type_file1 == 'prot' and type_file2 == 'prot' 
    factory_f1 = Bio::Blast.local('blastp', "./Databases/#{file1_db.to_s}")
    factory_f2 = Bio::Blast.local('blastp', "./Databases/#{file2_db.to_s}")

end
puts "\nRunning #{factory_f1.program} with file 1 and #{factory_f2.program} with file 2..."

# We create Bio::FastaFormat objects for each of the files
fasta_file1 = Bio::FastaFormat.open(file1)
fasta_file2 = Bio::FastaFormat.open(file2)

# We create a hash to store the file2 sequences and make them accesible by its id.
file2_hash = Hash.new 
fasta_file2.each do |seq_f2|
  file2_hash[(seq_f2.entry_id).to_s] = (seq_f2.seq).to_s 
end

# Count number orthologues found
sum = 0
evalue = 10**-6

puts "Searching for reciprocal best hits..."
fasta_file1.each do |seq_f1| # We iterate over each sequence in the file1 
  # We extract the ID of each sequence in file1
  file1_id = (seq_f1.entry_id).to_s 
  # Run the actual BLAST by querying the factory of file 2 with each of the sequences of the file1
  report_f2 = factory_f2.query(seq_f1)
  
  # Only keep the BLAST result if it passes the evalue filter
  if report_f2.hits[0] && report_f2.hits[0].evalue <= evalue  
    # Extract the ID of each sequence in file2
    file2_id = (report_f2.hits[0].definition.match(/(\w+\.\w+)|/)).to_s 
    # Run the actual BLAST by querying the factory of file 1 with each of the sequences of the file2
    # To get the sequence of each file2_id, we look at the hash we create previously where we linked each ID with its sequence
    report_f1 = factory_f1.query(">#{file2_id}\n#{file2_hash[file2_id]}")
    
    # Only keep the BLAST result if it passes the evalue filter
    if report_f1.hits[0] && report_f1.hits[0].evalue <= evalue  
      # Store the ID that will match with the ID in file1
      f1_match = (report_f1.hits[0].definition.match(/(\w+\.\w+)|/)).to_s 
      # If the match and the file1_file ID match = reciprocal best hit 
      if file1_id == f1_match
        #puts "\t#{file1_id}\t\t----------\t\t#{file2_id}" 
        output.puts "\t#{file1_id}\t\t\t----------\t\t#{file2_id}"
        sum += 1
        if sum%100==0
            puts "Orthologues found: #{sum}"
            puts "Resuming search..."
        end      
        
      end
    end         
  end
end

puts "\nDone"
puts "---------------------------------------------\n"
puts "TOTAL OF ORTHOLOGUES FOUND: #{sum}"
output.puts "\nTOTAL OF ORTHOLOGUES FOUND: #{sum}"
output.close

system("rm -r Databases")