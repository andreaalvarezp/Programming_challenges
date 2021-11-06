
# Assignment2: Gene.rb
# Author: Andrea Álvarez Pérez

# Modules

require 'json'
require 'net/http'
require './Interaction.rb'

class Gene
  
  attr_accessor :gene_id
  attr_accessor :prot_id
  attr_accessor :kegg_id          
  attr_accessor :kegg_path
  attr_accessor :go_id             
  attr_accessor :go_term
  
  @@total_genes = Hash.new      # Save all genes
    
  def initialize (params = {})
    
    @gene_id = params.fetch(:gene_id, "AT0G00000")
    @prot_id = params.fetch(:prot_id, "000000")
    
    # We do the annotations in arrays in case there are more than one
    
    @kegg_id = params.fetch(:kegg_id, Array.new)
    @kegg_path = params.fetch(:kegg_path, Array.new)
    @go_id = params.fetch(:go_id, Array.new)
    @go_term = params.fetch(:go_term, Array.new)
    
    # Save all genes to posterior annotation associated with its prot_id
    
    @@total_genes[prot_id] = self       
          
  end
  
  # Function to load input file
  
  def self.load_file(file, out)
    
    prueba = File.open(file, "r")
    
    prueba.each_line do |line|
      line.delete!("\n")
      puts;puts "AGI Locus: #{line} loaded"
      
      # Get prot_id for each AGI Locus code
      prot_id = self.get_prot_id(line)
      
      # Create a new Gene object
      
      Gene.new(
            :gene_id => line,
            :prot_id => prot_id
            )
      
      # Create Protein object and pass the gene_id information
      # Start in depth = 0
      Interaction.get_prot(prot_id, 0, gene_id = line)

    end

    prueba.close
    
  end
  
  # Function to call the method which stores all gene objects
  
  def self.genes
    
    return @@total_genes
    
  end
  
  # Function to get the prot_id 
  
  def self.get_prot_id(gene_id)
    
    response = fetch("http://togows.org/entry/ebi-uniprot/#{gene_id}/entry_id.json")
    
    if response
      
      data = JSON.parse(response)
      
      return data[0]
    else
      puts "Web call failed in function get_prot_id - see STDERR for details..."
    end
    
  end
  
  def annotate_data
    # Annotate GO Biological process term and KEGG pathway
    # I have to call the gene_id that is entering the class
    
    response_GO = fetch ("http://togows.org/entry/ebi-uniprot/#{self.gene_id}/dr.json")
    response_KEGG = fetch ("http://togows.org/entry/kegg-genes/ath:#{self.gene_id}/pathways.json")
    
    if response_GO || response_KEGG

      data_GO = JSON.parse(response_GO.body)
      data_KEGG = JSON.parse(response_KEGG.body)
      
      argument = "P:" # Biological process GO Terms only
  
      # GO data information is stored in an array
      data_GO[0]["GO"].each do |g|
        if (g[1].match("#{argument}"))    # Only biological process GO
          # Store the information in the corresponding arrays
          self.go_id << g[0]
          self.go_term << g[1]
        else
          next
        end
      end
        
      # KEGG data information is stored in a dictionary  
      data_KEGG[0].each do |kegg_id, kegg_path|
        # Store the information in the corresponding arrays
        self.kegg_id << kegg_id
        self.kegg_path << kegg_path
      end
    else
      puts "Web call failed in function annotate_data - see STDERR for details..."
    end
  end
    
end