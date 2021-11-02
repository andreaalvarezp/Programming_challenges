
# Assignment2: Gene.rb
# Author: Andrea Álvarez Pérez


require 'json'
require 'net/http'
require './Interaction.rb'

class Gene
  
  attr_accessor :gene_id
  attr_accessor :prot_id
  attr_accessor :kegg_id           # kegg id value
  attr_accessor :kegg_path
  attr_accessor :go_id             # go id value
  attr_accessor :go_term
  
  @@total_genes = Hash.new      # Save all genes
    
  def initialize (params = {})
    
    @gene_id = params.fetch(:gene_id, "AT0G00000")
    @prot_id = params.fetch(:prot_id, "000000")
    @kegg_id = params.fetch(:kegg_id, Array.new)
    @kegg_path = params.fetch(:kegg_path, Array.new)
    @go_id = params.fetch(:go_id, Array.new)
    @go_term = params.fetch(:go_term, Array.new)
    
    @@total_genes[prot_id] = self       # Save all genes to posterior annotation
          
  end
  
  def self.load_file(file, out)
    
    prueba = File.open(file, "r")
    
    prueba.each_line do |line|
      line.delete!("\n")
      puts;puts "AGI Locus: #{line} loaded"
      prot_id = self.get_prot_id(line)

      #File.open(out, "a"){ |f| f.write("\r\nAGI Locus: #{line}")}
      
      # Create a new Gene object
      
      Gene.new(
            :gene_id => line,
            :prot_id => prot_id
            )
      
      # Create Protein object and pass the gene_id information
      #puts "Get_prot object called by Gene class..."
      Interaction.get_prot(prot_id, 0, gene_id = line)
      #self.annotate_data(line, out)
    end

    prueba.close
    
  end
  
  def self.genes
    
    return @@total_genes
    
  end
  
  def self.get_prot_id(gene_id)
    
    response = fetch("http://togows.org/entry/ebi-uniprot/#{gene_id}/entry_id.json")
    
    if response
      
      data = JSON.parse(response)
      
      #puts "Prot ID of AGI Locus #{gene_id}: #{data[0]}"
      return data[0]
    else
      puts "Web call failed in function get_prot_id - see STDERR for details..."
    end
    
  end
  
  def annotate_data
    # Annotate GO Biological process term and KEGG pathway
    # I have to call the gene_id that is entering the class
    
    #puts "Annotating gene #{self.gene_id}..."
    response_GO = fetch ("http://togows.org/entry/ebi-uniprot/#{self.gene_id}/dr.json")
    response_KEGG = fetch ("http://togows.org/entry/kegg-genes/ath:#{self.gene_id}/pathways.json")
    
    if response_GO || response_KEGG

      data_GO = JSON.parse(response_GO.body)
      data_KEGG = JSON.parse(response_KEGG.body)
      
      argument = "P:" # Biological process GO Terms only
  
      # GO data information is stored in a list
      data_GO[0]["GO"].each do |g|
        #puts "Biological process GO terms:"
        if (g[1].match("#{argument}"))    
          self.go_id << g[0]
          self.go_term << g[1]
          #File.open(out, "a"){ |f| f.write("\r\n\tGO ID: #{go_id}\t\tGO TERM: #{go_term}")}
        else
          next
        end
      end
        
      # KEGG data information is stored in a dictionary  
      data_KEGG[0].each do |kegg_id, kegg_path|
        self.kegg_id << kegg_id
        self.kegg_path << kegg_path
        #File.open(out, "a"){ |f| f.write("\r\n\tKEGG ID: #{kegg_id}\t\tKEGG PATHWAY NAME: #{kegg_path}")}
      end
      #puts "Done!";puts
    else
      puts "Web call failed in function annotate_data - see STDERR for details..."
    end
  end
    
end