# => DNA SEQUENCE ALIGNMENT
#
# length of the first sequence, dna1, sets the number of rows
# length of the second sequence, dna2, sets the number of columns
#
# EXAMPLE>
# dna1 = 'GAG'
# dna2 = 'TCTCTC'
# 
# Starting matrix is:
#
#   _TCTCTC
#   ~~~~~~~
# _|0000000
# G|0000000
# A|0000000
# G|0000000
#

# libraries
    # require 'nmatrix' # using SciRuby NMatrix (gem install nmatrix)
    require 'matrix' # normal ruby Matrix class

    def pretty_matrix(dna1, dna2, mat) # seq1, seq2, scored_matrix
      # This function gives out the scoring matrix in a format 
      # well readable in the console

        puts 'The score matrix:'
        dna2.each_char {|d| print '    ' + d + ''}
        puts

          i = 0
          r = 0
          print dna1[r]
        mat.each do |number|
          if number >= 0          # za znaminko
            print ' '
          end
          if number >= 0 and number < 10    # za nulu
            print '  '
          end
          if number >= 10 and number < 100  # za nulu
            print ' '
          end
          if number > -10 and number < 0
            print '  '
          end
          if number > -100 and number <= -10
            print ' '
          end
          print number.to_s + " "
          i+= 1
          if i == mat.column_size
            print "\n"
            r += 1
            print dna1[r]
            i = 0
          end
        end
      return
    end

# User input (alignment params)
    puts "Input the MATCH scoring coefficient>"; STDOUT.flush
        match_coeff = gets.chomp.to_i

    puts "Input the MISMATCH scoring coefficient>"; STDOUT.flush
        mismatch_coeff = gets.chomp.to_i

    puts "Input the GAP scoring coefficient>"; STDOUT.flush
        gap_coeff = gets.chomp.to_i

# User input (sequences)
    correctdna1 = false
    while correctdna1==false
    
      puts "Input first DNA sequence>"; STDOUT.flush
        dna1 = gets.chomp.upcase
        if ( dna1 =~ /^[CAGTcagt]+$/ )
           correctdna1=true
       else
           puts 'Wrong characters detected, use only GATC.'
        end
    end

    correctdna2 = false
    while correctdna2==false
      puts "Input second DNA sequence>"; STDOUT.flush
        dna2 = gets.chomp.upcase
        if ( dna2 =~ /^[CAGTcagt]+$/ )
           correctdna2=true
       else
           puts 'Wrong characters detected, use only GATC.'
        end
    end
    
# Score the matrix
  # 1) Take the sequences, reverse them, add gap to the beginning (- is gap)
    dna1 = dna1.reverse
    dna2 = dna2.reverse
    dna1.insert(0, '-')
    dna2.insert(0, '-')
  # 2) Construct the empty matrix (using Matrix class)
    mat = Matrix.build(dna1.length, dna2.length) {|row, col| 0 }

  # 3) Fill in the first col and row
    for n in 0..dna2.length-1 # columns
      mat.send(:[]=,0,n,n*gap_coeff)
    end
    for n in 0..dna1.length-1 # rows
      mat.send(:[]=,n,0,n*gap_coeff)
    end

  # 4) Compute the scores
  # for each cell, look at top, left and top-left cell
    for row in 1..mat.row_count-1
      for col in 1..mat.column_count-1
  # a) top cell value, add gap_coeff to it
      topval = mat.[](row-1,col)
      topval = topval + gap_coeff 
  # b) left cell value, add gap_coeff to it
      leftval = mat.[](row,col-1)
      leftval = leftval + gap_coeff
  # c) diagonal value, just take it over
      diagval = mat.[](row-1,col-1)
  # d) check for match or mismatch, add/subtract coefficient
      top = dna2.[](col)
      left = dna1.[](row)

      if top == left
        topval = topval + match_coeff
        leftval = leftval + match_coeff
        diagval = diagval + match_coeff
      else
        topval = topval + mismatch_coeff
        leftval = leftval + mismatch_coeff
        diagval = diagval + mismatch_coeff
      end
  # e) select the highest score and write it
        maxval = [topval,leftval,diagval].max
        mat.send(:[]=,row,col,maxval)
      end
    end

  # f) show matrix
puts pretty_matrix(dna1,dna2,mat)
  # ----------------------------------------------------------------  

# Needleman-Wunsh (global) or Smith-Waterman (local?)


## TODO code for choosing


# NEEDLEMAN-WUNSH --- GLOBAL ALIGNMENT
# G) find the optimal path
# a) init + start in the lower right corner
    g_path = []
    g_aligned1 = []
    g_aligned2 = []
    row = dna1.length-1 # starting row
    col = dna2.length-1 # starting column
    g_path << mat.[](row,col)

# b) iterate over the path until top-left
    while row>0 or col>0
        lv = mat.[](row,col-1) # left value
        tv = mat.[](row-1,col) # top value
        dv = mat.[](row-1,col-1) # diagonal value

        if row == 0             # can go only to the left
          # puts 'goin left!!'
          g_path << lv
          g_aligned1 << dna2.[](col)
          g_aligned2 << '-' #dna1.[](row)
          col -= 1; # row remains the same
        elsif col == 0          # can go only up
          # puts 'goin up!!!'
          g_path << tv
          g_aligned1 << '-' #dna2.[](col)
          g_aligned2 << dna1.[](row)
          row -= 1; # col remains the same
        else
          if dv == [lv,tv,dv].max
            # puts 'goin diag from row:' + row.to_s + ', col:' + col.to_s
            g_path << dv
            g_aligned1 << dna2.[](col)
            g_aligned2 << dna1.[](row)
            row -= 1; col -= 1
          elsif tv == [lv,tv,dv].max
            # puts 'goin up from row:' + row.to_s + ', col:' + col.to_s
            g_path << tv
            g_aligned1 << '-' #dna2.[](col)
            g_aligned2 << dna1.[](row)
            row -= 1; # col remains the same
          elsif lv == [lv,tv,dv].max
            # puts 'goin left from row:' + row.to_s + ', col:' + col.to_s 
            g_path << lv
            g_aligned1 << dna2.[](col)
            g_aligned2 << '-' #dna1.[](row)
            col -= 1; # row remains the same
        end 
      end
    end

# finish
  puts
  puts 'Path (with upward preference):'
  print g_path; puts
  puts
  puts 'DNA alignment:'
  print g_aligned1; puts
  print g_aligned2; puts
  puts
  print 'Score: '; print score = g_path.inject(0, :+); puts

############## END OF GLOBAL ALIGNMENT #############################

########### Smith-Waterman // local alignment
# 1) set negative values to zeros
  for r in 0..mat.row_count-1
    for c in 0..mat.column_count-1
      if (mat.[](r,c)) < 0.to_i
        mat.send(:[]=,r,c,0)
      end
    end
  end
# 2) find the maximum of the matrix
  max = mat.max
  max_index = [] #array for row, col index of max

  for r in 0..mat.row_count-1
    for c in 0..mat.column_count-1
      if mat.[](r,c) == max
        max_index[0] = r
        max_index[1] = c
        breakthru = true
      end
      break if breakthru == true
    end
    break if breakthru == true
  end

# 3) We've got the init address, now go back until we reach 0 value
  l_path = []
  l_aligned1 = []
  l_aligned2 = []
  row = max_index[0]
  col = max_index[1]

    while row>0 or col>0

      if mat.[](row,col) == 0
        break
      end

        lv = mat.[](row,col-1) # left value
        tv = mat.[](row-1,col) # top value
        dv = mat.[](row-1,col-1) # diagonal value

        if row == 0             # can go only to the left
          # puts 'goin left!!'
          l_path << lv
          l_aligned1 << dna2.[](col)
          l_aligned2 << '-' #dna1.[](row)
          col -= 1; # row remains the same
        elsif col == 0          # can go only up
          # puts 'goin up!!!'
          l_path << tv
          l_aligned1 << '-' #dna2.[](col)
          l_aligned2 << dna1.[](row)
          row -= 1; # col remains the same
        else
          if dv == [lv,tv,dv].max
            # puts 'goin diag from row:' + row.to_s + ', col:' + col.to_s
            l_path << dv
            l_aligned1 << dna2.[](col)
            l_aligned2 << dna1.[](row)
            row -= 1; col -= 1
          elsif tv == [lv,tv,dv].max
            # puts 'goin up from row:' + row.to_s + ', col:' + col.to_s
            l_path << tv
            l_aligned1 << '-' #dna2.[](col)
            l_aligned2 << dna1.[](row)
            row -= 1; # col remains the same
          elsif lv == [lv,tv,dv].max
            # puts 'goin left from row:' + row.to_s + ', col:' + col.to_s 
            l_path << lv
            l_aligned1 << dna2.[](col)
            l_aligned2 << '-' #dna1.[](row)
            col -= 1; # row remains the same
        end 
      end
    end


  pretty_matrix(dna1,dna2,mat)
  puts
  puts 'Path (with upward preference):'
  print l_path; puts
  puts
  puts 'DNA alignment:'
  print l_aligned1; puts
  print l_aligned2; puts
  puts
  print 'Score: '; print score = l_path.inject(0, :+); puts

######## END OF LOCAL ALIGNMENT ##########################################