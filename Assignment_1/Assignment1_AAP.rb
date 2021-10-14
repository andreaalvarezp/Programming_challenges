### ASSIGNMENT 1: CREATING OBJECTS
## Author: Andrea Álvarez Pérez

## Usage: $ ruby Assignment1_AAP.rb
## Once the code is executed, a file "new_stock_file.tsv" is created or rewrited (if it alredy exists).

class SeedStock # superclass
  
  attr_accessor :seed_stockID
  attr_accessor :geneID
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_rem
  attr_accessor :example
    
  def initialize (params = {}) # default values

    @seed_stockID = params.fetch(:seed_stockID, "000000")
    @geneID = params.fetch(:geneID, "unknown")
    @last_planted = params.fetch(:last_planted, "unknown")
    @storage = params.fetch(:storage, "unknown")
    @grams_rem = params.fetch(:grams_rem, "0")
    
  end
  
  # write the header in output file. If the file exists, rewrite it
  File.write("new_stock_file.tsv", ["Seed_Stock", "Mutant_Gene_ID", "Last_Planted", "Storage", "Grams_Remaining"].join("\t"), mode: "w") 
   
  def seed_grams_remain
    newvalue = grams_rem.to_i - 7 # convert string into integer and plant 7 seeds
    if newvalue <= 0 then
      # set value into 0 to avoid having negative numbers
      newvalue = 0
      File.open("new_stock_file.tsv", "a") { |f| f.write("\r\n#{seed_stockID}\t#{geneID}\t#{last_planted}\t#{storage}\t#{newvalue}\t") }
      return "WARNING: we have run out of Seed Stock #{seed_stockID}"
    else
      # write results in new_stock_file
      File.open("new_stock_file.tsv", "a") { |f| f.write("\r\n#{seed_stockID}\t#{geneID}\t#{last_planted}\t#{storage}\t#{newvalue}\t") }
      return
    end
  end
  
  
  def calculate f2_wild = 0, f2_P1 = 0, f2_P2 = 0, f2_P1P2 = 0, gene_name = 0, gene2 = 0
    # First we calculate the expected frequencies
    # 1. Sum of all individuals
    sum = f2_wild.to_f + f2_P1.to_f + f2_P2.to_f + f2_P1P2.to_f
    # 2. 9:3:3:1 proportion
    f2_wild_E = (sum*9)/16
    f2_P1_E = (sum*3)/16
    f2_P2_E = (sum*3)/16
    f2_P1P2_E = sum/16
    
    #3. Chi-square test: (O-E)^2 / E
    chi2 = ( (((f2_wild.to_f - f2_wild_E)**2) / f2_wild_E) + (((f2_P1.to_f - f2_P1_E)**2) / f2_P1_E) + (((f2_P2.to_f - f2_P2_E)**2) / f2_P2_E) + (((f2_P1P2.to_f - f2_P1P2_E)**2) / f2_P1P2_E) )
    
    #4. Cutoff chi2 value: alfa = 0.05 and df(number of classes - 1) = 4 - 1 (Chi-square distribution table)
    cutoff = 7.815
    
    #5. Compare with hypothesis
    if chi2 > cutoff
      puts "Recording: #{gene_name} is genetically linked to #{gene2} with chisquare score #{chi2}"   
      puts; puts
      puts "Final report:"
      puts 
      puts "#{gene_name} is linked to #{gene2}" 
      puts "#{gene2} is linked to #{gene_name}"
    end
    
  end

end 


class HybridCross < SeedStock
  
  attr_accessor :parent1
  attr_accessor :parent2
  attr_accessor :f2_wild
  attr_accessor :f2_P1
  attr_accessor :f2_P2
  attr_accessor :f2_P1P2
  
  def initialize (params = {}) # default values
    
    super(params)
    @parent1 = params.fetch(:parent1, "unknown")
    @parent2 = params.fetch(:parent2, "unknown")
    @f2_wild = params.fetch(:f2_wild, "000")
    @f2_P1 = params.fetch(:f2_P1, "000")
    @f2_P2 = params.fetch(:f2_P2, "000")
    @f2_P1P2 = params.fetch(:f2_P1P2, "000")
       
  end
  
  def calculate (gene_name, gene2)
    
    super f2_wild, f2_P1, f2_P2, f2_P1P2, gene_name, gene2 #variable values to pass to the main class
    
  end
  
end

class Gene < HybridCross
  
  attr_accessor :geneID 
  attr_accessor :gene_name
  attr_accessor :mutant_phenotype
  
  def initialize (params = {}) # default values

    super(params)
    @geneID = params.fetch(:geneID, "00000")
    @gene_name = params.fetch(:gene_name, "unknown")
    @mutant_phenotype = params.fetch(:mutant_phenotype, "unknown")
    
  end
  
  def calculate 

    super gene_name
    
  end
 
end

# =======================================================================================================

# Open files and fill class data

require "csv"

# read CSV files and store the information in 3 variables
tsv_seed = CSV.read("seed_stock_data.tsv", col_sep: "\t")[1 .. -1] # delete the header
tsv_gene = CSV.read("gene_information.tsv", col_sep: "\t")[1 .. -1] # delete the header
tsv_cross = CSV.read("cross_data.tsv", col_sep: "\t")[1 .. -1] # delete the header

for i in [0,1,2,3,4]
  p2 = SeedStock.new(
    :seed_stockID => tsv_seed[i][0],      
    :geneID => tsv_seed[i][1],                   
    :last_planted => tsv_seed[i][2], 
    :storage => tsv_seed[i][3],
    :grams_rem => tsv_seed[i][4]
    )
  puts p2.seed_grams_remain # execute first exercise function
end
  
for i in [0,1,2,3,4]
  p3 = HybridCross.new(    
    :parent1 => tsv_cross[i][0],                   
    :parent2 => tsv_cross[i][1], 
    :f2_wild => tsv_cross[i][2],
    :f2_P1 => tsv_cross[i][3],
    :f2_P2 => tsv_cross[i][4],
    :f2_P1P2 => tsv_cross[i][5],
    )
  p4 = Gene.new(    
    :geneID => tsv_gene[i][0],                   
    :gene_name => tsv_gene[i][1], 
    :mutant_phenotype => tsv_gene[i][2],
    )
  # places in gene_name list advance one position in order to associate with parent2 list
  if i == 4
      gene2 = tsv_gene[0][1]
  else
    i += 1
    gene2 = tsv_gene[i][1]
    i -= 1
  end
  # execute first exercise function: the arguments are the gene names associated with parent1 and parent2
  puts p3.calculate(p4.gene_name, gene2) 
end
