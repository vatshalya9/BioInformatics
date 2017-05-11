%------ PDP Using Skiena's Algorithm -----
%PDP 

function PDP(A, iftrace)
A = sort(A);                         % make L sorted for easy implementation
s = size(A);                         % n is the length of L
val = max(A);                        % width is the maximum element of L
A = A(A ~= val);                     % delete (width, L)
X = [0 val];
PLACE(A, X, val, 0 , iftrace);
end

    function PLACE(A, X, val, level, iftrace)
        indent = '----';
        for i=1:level
            indent = strcat(indent, '----');
        end
        if (isempty(A))   
            sol = ' ';
            XX = sort(X);
                for i=1:length(X)
                    sol = strcat(sol, sprintf(' %d', XX(i)));
                end
            disp(fprintf('%s Solution found:%s', indent, sol));
        else
            y = max(A);                 % y is the maximum element in L
            deltay = sort(delta(y, X)); % dy is delta(y,X)
             z = (val - y);
             dz = sort(delta(z,X));     %dz is delta(val-y,  X)
         if anyEq(y, A)                 %delta(y,x) in L
                if iftrace == 1
                    fprintf('%s try y=%d\n', indent,y);
                end
                X = cat(y,X);               % add y to x
                A = sort(y,A);
                A = removeAll(deltay, A);    % remove all delta(y,x) in L
                PLACE(A, X, val, level+1 , iftrace);
                %X = remove(y,A)  ; % remove y from x
                %A = sort(X,A); % add delta(y,x) to L
            else
                if iftrace == 1
                    fprintf('%s try y = %d FAILS\n', indent, y);
                end
                
                if anyEq(dz, A)             %delta(y,x) in L
                    if iftrace == 1
                    fprintf('%s try z=%d\n', indent,dz);
                    end
                
                X = cat(z,X); 
                A = sort(z,A);
                A = removeAll(deltay, A);   % remove all delta(y,x) in L
                PLACE(A, X, val, level+1 , iftrace);
                %X = remove(Z,A)  ; % remove y from x
                %A = sort(Z,A); % add delta(y,x) to L
            else
                if iftrace == 1
                    fprintf('%s try z = %d FAILS\n', indent, dz);
                end
            end
         end
end
        function delta = delta(y, X)
            delta = abs(X-y);
        end
        function rem = removeAll(y, X)
            x = X;
            for i = 1:size(A)
              index = find(A(i) == X);
              if (size(index) >= 1)
                  X(index) = [];
              end
            end    
            rem = x;
        end
        function sub = anyEq(y, X)
        if (size(y) > size(X))
            sub = 0;
            return;
         end
            A = sort(X);
            B = sort(y);
            found = false;

            a = 1;
            b = 1;
            C = [];
            while a <= length(A) && b <= length(B) %Comparing the length of arrays with variable defined.
                    if (A(a) < B(b))
                        a = a + 1;                  % incrementing a value

                    else if (A(a) > B(b))
                            b = b + 1;              %incrementing b value

                        else 
                            C = [C A(a)]; %adding all the array elements to another array
                            a = a + 1;
                            b = b + 1;
                        end
                    end
            end


         if (isequal(C, A)) % checking whether equal or not
            sub = true;
         else
            sub = false;

         end

    end
    end