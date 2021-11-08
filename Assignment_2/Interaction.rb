
# Assignment2: Interaction_protein.rb
# Author: Andrea Álvarez Pérez

# Modules

require 'json'
require 'net/http'
require './Gene.rb'

class Interaction
    
    attr_accessor :prot_id        # UniProt ID
    attr_accessor :intact_id      # IntAct ID
    attr_accessor :network        # Network ID
    
    @@prot_intact = Hash.new      # Not all the protein objects are going to have an IntAct ID, we store here the ones with it
    @@interactions = false        # Accumulates the interactions by pairs of IntAct Codes
    
    def initialize(params = {})
      
      @prot_id = params.fetch(:prot_id, "XXXXXX") 
      @intact_id = params.fetch(:intact_id, nil)
      @network = params.fetch(:network, nil)

      # Only add new object of there is an intact ID present
      if intact_id
        @@prot_intact[intact_id] = self
      end
      
      if @@interactions == false
          @@interactions = Array.new        # Array containing all interactions, we want it to be acumulative
          # we only start this method once, and store all interactions of all genes
      end

      
    end
    
    # Function to return to main script all protein objects created to posterior network building
    
    def self.return_method
        
        return @@prot_intact
        
    end
    
    # Function to create protein object and store prot_id and intact_id (if exists) for each of them
    
    def self.get_prot(prot_id, depth, gene_id = nil, intact_id = nil)
      
      if not intact_id                                  # If the gene_id i given, I search for Intact_id
        intact_id = self.get_intact_id(gene_id)
      end
      
      if intact_id && (depth.to_f < $max_depth)         # If the IntAct is given, I search for interactions
        new_interaction(intact_id, depth)                               
      end

      Interaction.new(
            :prot_id => prot_id,
            :intact_id => intact_id,
            :network => nil 
            )
      
      depth += 1                                        # Go ahead with the network depth
 
    end
    
    # Function to look if a protein object alredy exist in order to avoid repetitions in networks.
    # We control if the protein has an intact_id associated
    
    def self.prot_exist(intact_id)

        if @@prot_intact.has_key?(intact_id)
          return true
        
        else
          return false
        
        end
          
    end
    
    # Function to get intact_id from the AGI locus code. I chose to use togows.
    
    def self.get_intact_id(gene_id)
        
        
        # Get intact_id from gene_id
        response = fetch("http://togows.org/entry/ebi-uniprot/#{gene_id}/dr.json")
        
        if response
  
          data = JSON.parse(response)
      
          if data[0]['IntAct']
            return data[0]['IntAct'][0][0]  # IntAct ID result
          else
            return nil                      # Empty result (not stored in @@prot_intact)
          end
        else
            puts "Web call failed in function get_intact_id - see STDERR for details..."
        end
      
    end
    
    # Function to convert intact_id to prot_id in order to fill protein objects attributes
    
    def self.convert_intact_protid(intact_id)
        
        #Get prot_id from intact_id
        response = fetch("http://togows.org/entry/ebi-uniprot/#{intact_id}/entry_id.json")
        
        if response

            data = JSON.parse(response)
            prot_id = data[0]
            
            return prot_id
        else
            puts "Web call failed in convert_intact_protid - see STDERR for details..."
        end
  
    end
    
    # Function called to create new protein objects according to interactions found.
    
    def self.new_interaction(intact_id, depth)
        
        # I have to create the proteins with the elements of the network I got in get_interactions
        interactions = get_interactions(intact_id)
        depth += 1
        
        if interactions                 # I get the rest of the interactions
            intact = intact_id          # I don't want my IntAct ID to change when I create a new protein object
            
            interactions.each do |prot1, prot2|
                # Check if the added protein is the first or the second in the array and if it doesn't exist as an object
                
                if prot1 == intact && (not prot_exist(prot2))
                    prot2_id = convert_intact_protid(prot2)
                    # Create new protein object
                    self.get_prot(prot2_id, depth, gene_id = nil, intact_id = prot2)    # search for interactions at 2º depth level
                    
                elsif prot2 == intact && (not prot_exist(prot1))    
                    prot1_id = convert_intact_protid(prot1)
                    # Create new protein object
                    self.get_prot(prot1_id, depth, gene_id = nil, intact_id = prot1)     # search for interactions at 2º depth level
                    
                end
            end

        end

    end
    
    # Function to get the interactions of a certain IntAct code. I chose to use EBI's REST API.
  
    def self.get_interactions(intact_id)
        
        # This function retrieves the IntAct URL with the interactions of my protein.
        
        puts "Fetching for interactions for IntAct ID #{intact_id}..."
        response = fetch ("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{intact_id}?format=tab25")
        
        if response

            data = response.to_s.split("\n")
            
            # create an array to save the protein interactions by pairs
            interactions = Array.new
          
            data.each do |line|
                          
                uni = line.split("\t")
               
                uniprot1 = uni[0].split(":")[1]                     # First protein 
                uniprot2 = uni[1].split(":")[1]                     # Second protein
                score = uni[14].sub(/intact-miscore:/, "").to_f     # Quality score
    
                ## FILTERS
                
                if score < $misscore.to_f                                                  # Quality score > 0.45 (middle-quality) or 0.3 (middle-low quality)
                    next
                elsif uniprot1 == uniprot2
                    next
                elsif uni[9] =~ /taxid:3702/ && uni[10] =~ /taxid:3702/         # Species: Arabidopsis
                    if uniprot1 < uniprot2                                      # Alphabetical order
                        if not @@interactions.include?([uniprot1, uniprot2])    # Avoid repetitions
                            @@interactions << [uniprot1, uniprot2]
                            interactions << [uniprot1, uniprot2]
                        end
                    else                                                        # No alphabetical order
                        if not @@interactions.include?([uniprot2, uniprot1])    # Avoid repetitions
                            @@interactions << [uniprot2, uniprot1]
                            interactions << [uniprot2, uniprot1]
                        end
                    end
                end
            end
            
            # List in global variables all the interactions and proteins with IntAct ID to use them in main script
            $INT = @@interactions
           
            return interactions        # proteins which belong to the same network
        else
            puts "Web call failed in function get_interactions - see STDERR for details..."
        end

    end
  
end


