## ----------------------
## ASSIGNMENT 3
## Andrea Álvarez Pérez
## ----------------------

## USAGE: $ ruby Assignment3.rb <input file>

#In GFF3 format, start coordinate should always be less or equal to end coordinate.
#When doing the "resta", number stay negative, so that why we insert them in the index [0].all?

#match_end = desde la posición de inicio de mi pattern, le sumo la longitud del pattern CTTCTT
#y como quiero saber la coordenada de la ultima letra, le resto 1.


# 1. sacar las secuencias de los genes de la lista

require 'rest-client'
require 'bio'

def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end 

def gene_seq(gene_id)

  response = fetch ("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id}")
  
  if response
    body = response.body
    if body =~ /ID/ 
  
      # save the sequence in a Bio::EMBL object
      entry = Bio::EMBL.new(body)
      # convert database entry to biosequence
      nucl_seq = entry.to_biosequence
      
      return nucl_seq
    else
      begin  # Handle web errors
        puts "Unable to connect to database [ensemblgenomesgene]"  
        raise "This is an error" 
      rescue 
        puts "Loading next AGI Locus code..." 
      end
    end
  else
    
    puts "The Web call failed: #{url} - see STDERR for details..."
  end

end

## search for sequence pattern and determine which matches are inside of the exon.

def exon_pattern_pos(exon_id, matches, len_seq, exon_position, strand, len_pattern, source)
  
  # search for the given pattern in the exons 
  # IMPORTANT: we create a Hash per strand
  pattern_in_exon = Hash.new
  
  if strand == '-'
    
    #puts "matches: #{matches}"
    matches.each do |match|
      # position of the last nucleotide of the pattern in the exon (length = number of hits)
      match_end = match + len_pattern - 1

      # Check if the pattern all is inside the exon
      if (match >= exon_position[0].to_i) && (match <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)

        # Convert the positions to work with the supposed reverse strand positions
        m_end = len_seq - match_end 
        m_start = len_seq - match
        # start coordinate should be always less or equal than end coordinate (we are un reverse strand)
        #puts "[m_end, m_start] = [#{m_end}, #{m_start}]"
        pattern_in_exon[[m_end, m_start]] = [source, exon_id, '-']
      end
        
     end 

  elsif strand == '+'
    
    matches.each do |match|
        
      match_end = match + len_pattern - 1 
      
      # Check if all the pattern is inside the exon sequence
      if (match >= exon_position[0].to_i) && (match <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
        # The condition is established to see whether the target is inside the exon
        pattern_in_exon[[match, match_end]] = [source, exon_id, '+']
      end
        
     end    
    
  end
  
  if pattern_in_exon
    return pattern_in_exon
  end
end

# Function to search the pattern in all sequences of the gene list
# Then this sequences are filtered whether is an exon or not.
# For each exon with a pattern match we take: [start, end] position of the exon in the se

def find_pattern(nucl_seq, pattern, len_pattern, source)
    
  # Length of the nucleotide sequence
  len_seq = nucl_seq.length() 
  
  # Hash with all the positions os all the sequences found in the input gene list
  pattern_ex_positions = Hash.new   
  
  # 1. Make arrays with initial position of the pattern in the foward and reverse strands (length = number of matches)
  # Regexp.last_match.begin tells me the starting and ending coordinates of the match in the string
  foward_matches = nucl_seq.gsub(/#{pattern}/).map{Regexp.last_match.begin(0)}
  # reverse_matches: start position of the pattern in the reverse complement strand 
  reverse_matches = nucl_seq.reverse_complement.gsub(/#{pattern}/).map{Regexp.last_match.begin(0)}
  
  # 2. Take all the exons and their positions
  nucl_seq.features.each do |feature|
    
    if feature.feature == "gene" 
      feature.qualifiers.each do |qual|
        if qual.qualifier == "note"
          source << qual.value.split(":")[1].split(";")[0].delete(" ") 
        end          
      end
    end
    
    source.each do |s|
    # .feature option says if we have an exon
      if feature.feature == "exon" 

        # we need to know the strand where the exon is coded [start, end]
        position = feature.position

        # exon_id is stored in the VALUE qualifier of feature.
        exon_id = feature.qualifiers[0].value.gsub('exon_id=', '')
        
        case
          when position =~ /[A-Z]/ ##### QUITO O NO LOS QUE TIENEN LETRAS MAYUSCULAS
            next
          when position !~ /[A-Z]/
        
            # if the exon is in the reverse strand:
            if position =~ /complement/ 
              #puts "Exon in reverse strand: "
    
              # delete "complement" and separate positions by .. (start, end)
              # if don't delete position with letters: position = position.tr('complement()', '').split(':')[1].split('..')
              position = position.tr('complement()', '').split('..')
              position_reverse = []
    
              # coordinates of the exon in complement strand are not coded properly
              position.each do |pos|
                # full sequence length - position in the reverse strand
                position_reverse.insert(0, len_seq - pos.to_i) 
                # repeat for start and end positions of the exon
              end
              
              # 3. Compare if pattern matches appear to be in any of the exons
              target_pos_in_exon = exon_pattern_pos(exon_id, reverse_matches, len_seq, position_reverse, '-', len_pattern, s)
    
              # 4. If the patterns found are inside the exons, we add them to our hash
              if target_pos_in_exon 
               pattern_ex_positions = pattern_ex_positions.merge(target_pos_in_exon)
              end
    
            # exon in the foward strand
            else
    
              # separate positions by .. (start, end)
              position_foward = position.split('..') 
              
              # 3. Compare if pattern matches appear to be in any of the exons
              target_pos_in_exon = exon_pattern_pos(exon_id, foward_matches, len_seq, position_foward, '+', len_pattern, s)
    
               # 4. If the patterns found are inside the exons, we add them to our hash
              if target_pos_in_exon 
               pattern_ex_positions = pattern_ex_positions.merge(target_pos_in_exon)
              end
            end
        end
      else
        next
      end
    end
  end
  # pattern_ex_positions hash structure: [start, end] => [exon_id, strand]
  return pattern_ex_positions
  
end

# 4a. Write Features into a GFF3 file
def write_file(gene_id, exon_features, outfile1)
  
  File.open(outfile1, "a") do |gff|
  
    exon_features.each do |feat|
      # Start and end positions
      pos_start = feat.position.split("..")[0]
      pos_end = feat.position.split("..")[1]

      # I save the values of each qulifier in an array to write in the output file
      values = Array.new
      feat.qualifiers.each do |qual|
        values <<  qual.value
      end
      
      # Column 3 (type): based on existing SOFA sequence ontology, I chose 'exon_region' with SOFA ID SO:0000852
      # ref: http://www.sequenceontology.org/so_wiki/index.php/Category:SO:0000852_%21_exon_region
      gff.puts "#{gene_id}\t#{values[0]}\t#{feat.feature}\t#{pos_start}\t#{pos_end}\t.\t#{values[2]}\t.\tID=#{values[1]};Note=pattern_#{values[3]}"

    end  
  end

end

# Add features to the EnsEMBL Sequence object and loop over each one of your CTTCTT features (using the #features method)
def add_features(gene_id, pattern_hits, nucl_seq, pattern, outfile1)
  
  exon_features = Array.new
  
  # Separate hash key::value
  pattern_hits.each do |pos, prop| 
    
    #Bio::Feature: feature, position, qualifiers
    # Output structure: http://bioruby.org/rdoc/Bio/Features.html
    
    # Check with feat.feature
    new_feat = Bio::Feature.new("SO:0000852", "#{pos[0]}..#{pos[1]}")
    
    # Check with feat.qualifiers
    new_feat.append(Bio::Feature::Qualifier.new('source', "#{prop[0]}"))
    new_feat.append(Bio::Feature::Qualifier.new('exon_id', "#{prop[1]}"))
    new_feat.append(Bio::Feature::Qualifier.new('strand', "#{prop[2]}"))
    new_feat.append(Bio::Feature::Qualifier.new('repeat_motif', "#{pattern.upcase}"))

    exon_features << new_feat
    #$FEATURES = exon_features
  end
  
  puts exon_features
  
  write_file(gene_id, exon_features, outfile1)
  
end


input = ARGV

if input.length != 2 
  puts "\r\n------------------------------------------------------------------------------------------"
  puts "\tERROR: Incorrect number of arguments"
  puts "\tUSAGE: $ruby Assignment3.rb <input file name> <output file name>"
  puts "------------------------------------------------------------------------------------------";puts
  exit 1
end

pattern = "CTTCTT"
len_pattern = pattern.length()

outfile1 = input[1]

if File.exist?(outfile1)
  File.delete(outfile1)
end

# Files
File.open(outfile1, "w") do |gff|
  gff.puts "##gff-version 3"
end

puts "\n---------------------------------------------"
puts "Assignment 3"
puts "GFF feature files and visualization"
puts "Andrea Álvarez Pérez"
puts "---------------------------------------------\n"

inputfile = input[0]
genes = Array.new
puts "Loading gene list..."
File.readlines(inputfile).each do |line|
  genes << line.delete("\n")
end

no_pattern_genes = Array.new
source = Array.new

puts "Searching for sequence and features..."
puts "Writing GFF3 output file..."
genes.each do |gene_id|
  
  nucl_seq = gene_seq(gene_id)

  if nucl_seq

    pattern_hits = find_pattern(nucl_seq, pattern.downcase, len_pattern, source)
    #puts pattern_hits
    
    # If the pattern_ex_positions is empty, we add the gene to a list
    if pattern_hits.empty?
      no_pattern_genes << gene_id
      
    # if the gene has the pattern, we search for the features
    else  
      add_features(gene_id, pattern_hits, nucl_seq, pattern, outfile1) # We create new features and add them to each seq_obj
    
    end  
    
  else
    next
  end
end
puts "Done!"
#4.b Output a report showing which genes on your list do NOT have exons with the CTTCTT repeat
puts "\n-----------------------------------------------------------------------------------------------------"
puts "Genes on the list '#{inputfile}' which do NOT have exons with the #{pattern.upcase} repeat:"
no_pattern_genes.each do |gene_id|
  if gene_id.empty?
    next
  else
    puts "- Gene ID: #{gene_id}"
  end
end
puts "-----------------------------------------------------------------------------------------------------"
puts "\nResults loaded in files:"
puts "- #{outfile1}";puts

  