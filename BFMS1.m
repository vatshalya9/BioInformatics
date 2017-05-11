%---------------BRUTE FORCE MOTIF PROBLEM--------------%
function BFMS1(DNA, l)

  timeStart = tic;
  [~, bestScore, bestMotif, bestWord, ~] = BruteForceMotifSearch(DNA, l);  
  
  fprintf('Update: Best Score = %d at Best Motif = %s and Best Word = %s\n', bestScore, mat2str(bestMotif), bestWord);
  timeElapsed = toc(timeStart);
  fprintf('Time taken for BFMS1 is : %f seconds\n', timeElapsed);
end


% This is a function to find Best score for the Motif finding problem.

function [cstr, bestScore, bestMotif, bestWord, best] = BruteForceMotifSearch(DNA, l)

  t = size(DNA, 1); % gives number of rows in DNA.
  n = length(DNA); % gives the number of columns in the DNA 
  v = n - l + 1; 
  i = 1;
  s = [];
  best = [];
  
  nucleotides = getAlphabets(t, DNA);
  disp(nucleotides);
  % initializing values 
  for j = 1 : t
    s = [s i];
  end
  init = s;  
  
  [score, word] = scoreFunc(s, l, nucleotides, DNA);
  cstr = word;
  fprintf('Consensus string = %s at s = %s\n', mat2str(s), cstr);
  bestScore = score;
  bestMotif = s;
  bestWord = word;
  while true  

      s = nextLeaf(s, t, v);    

      if (isequal(s, init))
          return;
      end
      
      % find the score.
      [score, word] = scoreFunc(s, l, nucleotides, DNA);
      fprintf('Consensus string = %s at s = %s\n', mat2str(s), cstr);
      if (score > bestScore)      
          bestScore = score;
          bestMotif = s;
          bestWord = word;
      end    
   
  end % end of while

end % end of function

% function to find the next search tree motif.
function next = nextLeaf(a, L, k)
    next = [];
    for i = L : -1 : 1     
        if a(i) < k
            a(i) = a(i) + 1;
            next = a;
            return;
        else
            a(i) = 1;
            next = a;
        end       
    end
end

% function to find the score of DNA for a value of Motif
function [score, word] = scoreFunc(s, l, nucleotides, DNA)

    alignment = createAlignMatrix(s, DNA, l);
    profile = createProfile(alignment, nucleotides);
    word = '';    
    [m, index] = max(profile);
    score = sum(m,2);
        
    for i = 1 : length(index)
        word = strcat(word, nucleotides(index(i))); 
    end
    
end


function alignment = createAlignMatrix(s, DNA, l)
    alignment = [];
    for i = 1 : length(s)          
        str = DNA(i,:);
        substr = str(s(i):s(i)+l-1);
        alignment(i, :) = substr;
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

function alphabet = getAlphabets(t, DNA)
    alphabet = [];
    for i = 1 : t
        temp = unique(DNA(i, :));
        if (length(alphabet) < length(temp))
            alphabet = temp;
        end
    end
end