#
# => DNA SEQUENCE ALIGNMENT
#
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
#
#


# libraries
    # require 'nmatrix' # using SciRuby NMatrix (gem install nmatrix)
    require 'matrix' # normal ruby Matrix class

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
    
# Needleman-Wunsch algorithm

  # 1) Take the sequences, add gap to the beginning (- is gap)
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

  # jedeme po radcich
      actual_cell_val = mat.[](row,col)

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


  # f) THIS WHOLE SHIT ONLY SHOWS PRETTY MATRIX ------------------
    puts 'The score matrix:'
    # puts mat.to_a.map(&:inspect) # shows the score matrix

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
  # ----------------------------------------------------------------  


# 5) find the optimal path
# 
# a) init + start in the lower right corner
    path = []
    aligned1 = []
    aligned2 = []
    row = dna1.length-1 # starting row
    col = dna2.length-1 # starting column
    path << mat.[](row,col)

# b) iterate over the path until top-left

    while row>0 or col>0


        lv = mat.[](row,col-1) # left value
        tv = mat.[](row-1,col) # top value
        dv = mat.[](row-1,col-1) # diagonal value

        if row == 0             # can go only to the left
          # puts 'goin left!!'
          path << lv
          aligned1 << dna2.[](col)
          aligned2 << '-' #dna1.[](row)
          col -= 1; # row remains the same
        elsif col == 0          # can go only up
          # puts 'goin up!!!'
          path << tv
          aligned1 << '-' #dna2.[](col)
          aligned2 << dna1.[](row)
          row -= 1; # col remains the same
        else
          if dv == [lv,tv,dv].max
            # puts 'goin diag from row:' + row.to_s + ', col:' + col.to_s
            path << dv
            aligned1 << dna2.[](col)
            aligned2 << dna1.[](row)
            row -= 1; col -= 1
          elsif tv == [lv,tv,dv].max
            # puts 'goin up from row:' + row.to_s + ', col:' + col.to_s
            path << tv
            aligned1 << '-' #dna2.[](col)
            aligned2 << dna1.[](row)
            row -= 1; # col remains the same
          elsif lv == [lv,tv,dv].max
            # puts 'goin left from row:' + row.to_s + ', col:' + col.to_s 
            path << lv
            aligned1 << dna2.[](col)
            aligned2 << '-' #dna1.[](row)
            col -= 1; # row remains the same
        end 
      end
    end



# c) finish
  puts
  puts 'Path (with upward preference):'
  print path; puts
  puts
  puts 'DNA alignment:'
  print aligned1.reverse; puts
  print aligned2.reverse; puts
  puts
  print 'Score: '; print score = path.inject(0, :+); puts
