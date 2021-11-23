## ----------------------
## ASSIGNMENT 3
## Andrea Álvarez Pérez
## ----------------------

## ------ FUNCTIONS ------

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
    # save the sequence in a Bio::EMBL object
    entry = Bio::EMBL.new(body)
    # convert database entry to biosequence
    nucl_seq = entry.to_biosequence
    return nucl_seq
  else
      puts "The Web call failed - see STDERR for details..." 
      return
  end
end

## Search for sequence pattern and determine which matches are inside of the exon.
def exon_pattern_pos(exon_id, matches, len_seq, exon_position, strand, len_pattern, source)
  
  # search for the given pattern in the exons. IMPORTANT: we create a Hash per strand
  pattern_in_exon = Hash.new
  
  # We add 3 atrributes to pattern_in_exon: source, exon_id, and strand direction
  if strand == '-'
    matches.each do |match|
      # position of the last nucleotide of the pattern in the exon (length = number of hits)
      # Example: if matches = [1233], match_end = 1238, position of the end of the pattern in the reverse strand
      match_end = match + len_pattern - 1
      # Check if the pattern all is inside the exon
      if (match >= exon_position[0].to_i) && (match <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
        # According to ensembl, positions columns are defined "with sequence numbering starting at 1", so if my pattern positions in
        # the reverse strand are [1233..1238], to write it in the file I have to convert to the corresponding position according to the foward strand
        m_end = len_seq - match_end 
        m_start = len_seq - match
        # start coordinate should be always less or equal than end coordinate (sideways)
        pattern_in_exon[[m_end, m_start]] = [source, exon_id, '-']
      end 
     end 
  elsif strand == '+'
    matches.each do |match|  
      match_end = match + len_pattern - 1 
      # Check if all the pattern is inside the exon sequence
      if (match >= exon_position[0].to_i) && (match <= exon_position[1].to_i) && (match_end >= exon_position[0].to_i) && (match_end <= exon_position[1].to_i)
        pattern_in_exon[[match, match_end]] = [source, exon_id, '+']
      end  
     end    
  end
  
  if pattern_in_exon
    return pattern_in_exon
  end
end

# Function to search the pattern in all sequences of the gene list. Then this sequences are filtered whether is an exon or not.
# For each exon with a pattern match we take: [start, end] position of the exon in the se
def find_pattern(nucl_seq, pattern, len_pattern, source)
    
  # Length of the nucleotide sequence
  len_seq = nucl_seq.length() 
  # Hash with all the positions os all the sequences found in the input gene list
  pattern_ex_positions = Hash.new   
  
  # 1. Make arrays with initial position of the pattern in the foward and reverse strands (length = number of matches)
  # Regexp.last_match.begin tells me the start coordinate of the match in the string
  foward_matches = nucl_seq.gsub(/#{pattern}/).map{Regexp.last_match.begin(0)}
  # reverse_matches: start position of the pattern in the reverse complement strand
  # Example: if reverse_matches = [1233], you can find repetition AAGAAG in position len_seq - 1233 - 6 (len_pattern) of the sequence given in the entry
  # which corresponds to position 1233 on the reverse strand
  reverse_matches = nucl_seq.reverse_complement.gsub(/#{pattern}/).map{Regexp.last_match.begin(0)}
  
  # 2. Take all the exons and their positions
  nucl_seq.features.each do |feature|
    # Determine source column
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
          when position =~ /[A-Z]/
            next
          when position !~ /[A-Z]/
            # if the exon is in the reverse strand:
            if position =~ /complement/  
              # delete "complement" and separate positions by .. (start, end)
              position = position.tr('complement()', '').split('..')
              position_reverse = []
              # coordinates of the exon in complement strand 
              position.each do |pos|
                # full sequence length - position in the reverse strand
                position_reverse.insert(0, len_seq - pos.to_i) 
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
  # pattern_ex_positions hash structure: [start, end] => [source, exon_id, strand]
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
      gff.puts "#{gene_id}\t#{values[0]}\t#{feat.feature}\t#{pos_start}\t#{pos_end}\t.\t#{values[2]}\t.\tID=#{values[1]};Note=pattern_#{values[3]}"
    end  
  end

end

# Add features to the EnsEMBL Sequence object and loop over each one of your CTTCTT features (using the #features method)
def add_features(gene_id, pattern_hits, nucl_seq, pattern, outfile1)
  
  exon_features = Array.new
  # Separate hash key::value
  pattern_hits.each do |pos, prop|    
    # Feature type: based on existing SOFA sequence ontology, I chose 'exon_region' with SOFA ID SO:0000852
    # ref: http://www.sequenceontology.org/so_wiki/index.php/Category:SO:0000852_%21_exon_region
    new_feat = Bio::Feature.new("SO:0000852", "#{pos[0]}..#{pos[1]}")
    # Check with feat.qualifiers
    new_feat.append(Bio::Feature::Qualifier.new('source', "#{prop[0]}"))
    new_feat.append(Bio::Feature::Qualifier.new('exon_id', "#{prop[1]}"))
    new_feat.append(Bio::Feature::Qualifier.new('strand', "#{prop[2]}"))
    new_feat.append(Bio::Feature::Qualifier.new('repeat_motif', "#{pattern.upcase}"))

    exon_features << new_feat
  end
  write_file(gene_id, exon_features, outfile1)
  
end

# Function that searches the nucl_seq coordinates inside the whole chromosome
def chromosome_coord(gene_id, nucl_seq, outfile)

  pr_acc = nucl_seq.primary_accession  
  if pr_acc
    File.open(outfile, "a") do |gff_chrom|
      # We are interested in changing columns 1, 4 and 5
      chrom = pr_acc.split(":")
      # Positions of the whole gene inside the chromosome
      gff_chrom.puts "#{chrom[2]}\t.\tgene\t#{chrom[3]}\t#{chrom[4]}\t.\t+\t.\tID=#{gene_id}"
      chr_num = chrom[2]
      pos1 = chrom[3]
      pos2 = chrom[4]
      return chr_num, pos1, pos2
    end
  else
    return false
  end
  
end

def pattern_in_chr(gene_id, pattern_hits, chr, outfile)
  
  File.open(outfile, "a") do |gff_chrom|
    # Structure: positions = [start, end]; exon_strand = [source, exon_id, strand]
    # We have alredy seen that the patterns are inside exons, so it's not neccessary to locate exon coordinates inside the chromosome
    pattern_hits.each do |positions, features|
      start_chr = chr[1].to_i + positions[0].to_i
      end_chr = chr[1].to_i + positions[1].to_i
      
      gff_chrom.puts "#{chr[0]}\t.\tnucleotide_motif\t#{start_chr}\t#{end_chr}\t.\t#{features[2]}\t.\tID=#{features[1]};Parent=#{gene_id}"
    end
  end
end

## ------ CODE ------

input = ARGV

if input.length != 1 
  puts "\r\n------------------------------------------------------------------------------------------"
  puts "\tERROR: Incorrect number of arguments"
  puts "\tUSAGE: $ruby Assignment3.rb <input file name> <output file name>"
  puts "------------------------------------------------------------------------------------------";puts
  exit 1
end

pattern = "CTTCTT"
len_pattern = pattern.length()

outfile1 = "gff3_genes.gff"
outfile2 = "gff3_chrom.gff"

# Files
File.open(outfile1, "w") do |gff|
  gff.puts "##gff-version 3"
end

File.open(outfile2, "w") do |gff_chrom|
  gff_chrom.puts "##gff-version 3"
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
    # If the pattern_ex_positions is empty, we add the gene to a list
    if pattern_hits.empty?
      no_pattern_genes << gene_id
    # if the gene has the pattern, we search for the features
    else  
      add_features(gene_id, pattern_hits, nucl_seq, pattern, outfile1) # We create new features and add them to each seq_obj
      chromosome = chromosome_coord(gene_id, nucl_seq, outfile2)
      pattern_in_chr(gene_id, pattern_hits, chromosome, outfile2)
    end  
  else
    next
  end
end

puts "Done!"
#4.b Output a report showing which genes on your list do NOT have exons with the CTTCTT repeat
puts "\nOUTPUT REPORT"
puts "-----------------------------------------------------------------------------------"
puts "Number of genes loaded: #{genes.length}"
puts "Genes on the list '#{inputfile}' which do NOT have exons with the #{pattern.upcase} repeat: #{no_pattern_genes.length}"

no_pattern_genes.each do |gene_id|
  if gene_id.empty?
    next
  else
    puts "- Gene ID: #{gene_id}"
  end
end
puts "-----------------------------------------------------------------------------------"
puts "Results loaded in files:"
puts "- #{outfile1}"
puts "- #{outfile2}"
puts "-----------------------------------------------------------------------------------";puts

  