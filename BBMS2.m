%--------------------------BRAND AND BOUND MEDIAN SEARCH-------------------------%
function BBMS2(DNA, l)

  timeStart = tic;
  [~, bestScore, bestMotif, bestWord, ~] = BranchAndBoundMotifSearch(DNA, l);  
  
  fprintf('Update: Best Score = %d at Best Motif = %s and Best Word = %s\n', bestScore, mat2str(bestMotif), bestWord);
  timeElapsed = toc(timeStart);
  fprintf('Time taken for BBMS2 is : %f seconds\n', timeElapsed);
end


function [cstring, bestScore, bestMotif, bestWord, best] = BranchAndBoundMotifSearch(DNA, l)

  t = size(DNA, 1); % gives number of rows in DNA.
  n = length(DNA); % gives the number of columns in the DNA 
  last = n - l + 1; 
  i = 1;
  s = [];
  best = [];
  bestScore = 0;
  
  
  nucleotides = getAlphabets(t, DNA);
  
  % initialize values for motif vector s = (1,1,....1)
  for j = 1 : t
    s = [s i];
  end 

  i = 1;
  %fprintf('t is : %d \n', t);
  while (i > 0)
     timeStart = tic;
     if (i < t) 
         
         [score, word] = scoreF(s, l, i, nucleotides, DNA);
         optimalscore = score + (t - i) * l;
         cstring = word;
         if (optimalscore < bestScore)
             [s,i] = byPass(s, i, t, last);
         else
             [s,i] = nextVertex(s, i, t, last); 
         end
         
     else       
         [score, word] = scoreF(s, l, t, nucleotides, DNA);
         cstring = word;
         if (score > bestScore)
              bestScore = score;
              bestMotif = s;
              bestWord = word;         
         end
         [s, i] = nextVertex(s, i, t, last);
     end
     timeElapsed = toc(timeStart);
     %fprintf('Time taken for each BBMS2 is : %f seconds\n', tElapsed);
  end

end

% function to find the next search tree motif.
function [a, val] = nextVertex(a, i, L, k)
    
    if (i < L)
        a(i+1) = 1;
        val = i + 1;
        return;
    else        
         for j = L : -1 : 1
            if (a(j) < k)
                a(j) = a(j) + 1;
                val = j;
                return;
            end
         end
    end
    val = 0;
    return;
    
end

function [a, val] = byPass(a, i, L, k)

    for j = i : -1 : 1
        if (a(j) < k)
            a(j) = a(j) + 1;
            val = j;
            return;
        end
    end
    val = 0;
    return;
    
end

% function to find the score of DNA for a value of Motif
function [score, word] = scoreF(s, l, r, nucleotides, DNA)
    alignment = createAlignMatrix(s, l, r, DNA);
    profile = createProfile(alignment, nucleotides);
        
    [m, index] = max(profile);
    score = sum(m,2);
    
    word = '';
    for i = 1 : length(index)
        word = strcat(word, nucleotides(index(i))); 
    end
    
end


function alignment = createAlignMatrix(s, l, r, DNA)
    alignment = [];
   
     for j = 1 : r       
        str = DNA(j, :);
        k = s(j);
        substr = str(k : k+l-1);
        alignment(j, :) = substr;

     end

end

function profile = createProfile(align, nucleotides)
    profile = [];  
    for i = 1 : length(nucleotides)
        for j = 1 : length(align(1, :))  
            profile(i, j) = sum(align(:, j) == nucleotides(i));
        end
    end
end

function alphabets = getAlphabets(t, DNA)
    alphabets = [];
    for i = 1 : t
        tmp = unique(DNA(i, :));
        if (length(alphabets) < length(tmp))
            alphabets = tmp;
        end
    end
end