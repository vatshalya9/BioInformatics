%---------------------------------Program for BRANCH AND BOUND MEDIAN SEARCH ---------------------%%
% Function BBMS
function BBMS4(DNA, l)

  timeStart = tic;
  [bestDistance, bestWord] = BranchAndBoundMedianSearch(DNA, l);  
  
  fprintf('Update: Best Distance = %d for Best Word = %s\n', bestDistance, bestWord);
  timeElapsed = toc(timeStart);
  fprintf('Time taken for BBMS4 is : %f seconds\n', timeElapsed);
end


function [bestDistance, bestWord] = BranchAndBoundMedianSearch(DNA, l)

    t = size(DNA, 1); % gives number of rows in DNA.
    
    k = 1;
    s = [];
    bestDistance = Inf;
    bestWord = '';
    
    nucleotides = getAlphabets(t, DNA);
    
    % initializing values 
    for j = 1 : l
        s = [s k];
    end

    i = 1;
    
    while (i > 0)
        
        if (i < l)
           
           prefix = '';
           for p = 1 : i
                q = s(p);
                prefix = strcat(prefix, nucleotides(q)); 
            end
            
            optimalDist = totalDistance(prefix, DNA, i);
            if (optimalDist > bestDistance)
               [s, i] = byPass(s, i, l, 4);
            else
                [s,i] = nextVertex(s, i, l, 4);
            end          
            
        else
            
            word = '';
            for p = 1 : l
                q = s(p);
                word = strcat(word, nucleotides(q)); 
            end
            
            totalDist = totalDistance(word, DNA, l);
            if (totalDist < bestDistance)
                bestDistance = totalDist;
                bestWord = word;
            end
            [s,i] = nextVertex(s, i, l, 4);        
            
        end
        
    end

end

% function for finding the next search tree motif.
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

% function for finding the score of DNA for a value of Motif
function dist = totalDistance(word, DNA, l)

    t = size(DNA, 1);
    n = length(DNA);
    min = Inf;
    hDist = 0;
    
    for j = 1 : t
        
        str = DNA(j, :);
     
        for k = 1 : n-l+1
            substr = str(k : k+l-1);
            diff = sum (substr ~= word);
            if (diff < min)
                min = diff;
            end
        end
        
        hDist = hDist + min;
    end
    dist = hDist;
    
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