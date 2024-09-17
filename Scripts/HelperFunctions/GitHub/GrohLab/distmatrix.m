function [ dist ] = distmatrix( P1,P2,n )
%DISTMATRIX Distance matrix between each point of vectors P1 and P2.
%   n is the norm from 1 to a reasonable number (not inf included).
if nargin == 2
    n = 2;
end
if nargin > 3
    error('This function just accepts 2 input vectors, and the norm order.')
elseif isnumeric(n)
    [points,dimension]=size(P1);
    [points2,dimension2]=size(P2);
    if dimension == dimension2
        dist = zeros(points,points2);
        for pf=1:points
            for pm=1:points2
                dist(pf,pm) = sum((P2(pm,:)-P1(pf,:)).^n)^(1/n);
            end
        end
    else
        error(['Vector 2 has ',num2str(dimension2),...
            ' dimensions, and vector 1 has ',num2str(dimension),...
            ' dimensions'])
    end
elseif ischar(n)
    switch n
        case 'pdf'
            % We don't mess with length of the pdf.
            [pdfs1,len1] = size(P1);
            [pdfs2,len2] = size(P2);
            if len1 == len2
                dist = zeros(pdfs1, pdfs2);
                for cpdf = 1:pdfs1
                    for cmpdf = 1:pdfs2
                        dist(cpdf,cmpdf) =...
                            KullbackLeiblerDivergence(P1(cpdf,:),P2(cmpdf,:));
                    end
                end
            else
                fprintf('Make sure that the distributions have equal length.\n');
                dist = NaN;
                return;
            end
        otherwise
            fprintf('''pdf'' for comparing probability density functions.\n')
            dist = NaN;
            return;
    end
end
end

