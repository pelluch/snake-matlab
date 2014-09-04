function Y = MonotoneConvHull(X)
%    Input: an iterable sequence of (x, y) pairs representing the points.
%   Output: a list of vertices of the convex hull in counter-clockwise order,
%      starting from the vertex with the lexicographically smallest coordinates.
%    Implements Andrew's monotone chain algorithm. O(n log n) complexity.


% Sort the points lexicographically (tuples are compared lexicographically).
% Remove duplicates to detect the case we have just one unique point.
% points = sort(unique(X,'rows'),1);
points = sortrows(X,1);
% Boring case: no points or a single point, possibly repeated multiple times.
if size(points,1) <= 1
    Y = points;
end

% Build lower hull
lower = [];
for p=1:size(points,1)
    if(size(lower,1) >=2)
        cp = crosspod(lower(end-1,:), lower(end,:), points(p,:));
    end
    while (size(lower,1) >= 2 && crosspod(lower(end-1,:), lower(end,:), points(p,:)) <= 0)
        lower = lower(1:end-1,:);
    end
    lower(end+1,:) = points(p,:);
end

% Build upper hull
upper = [];
for p=size(points,1):-1:1
    if(size(upper,1) >=2)
        cp = crosspod(upper(end-1,:), upper(end,:), points(p,:));
    end
    while (size(upper,1) >= 2 && crosspod(upper(end-1,:), upper(end,:), points(p,:)) <= 0)
        upper = upper(1:end-1,:);
    end
    upper(end+1,:) = points(p,:);
end

%Remove last point of each list ( its the same as the first of the other)
lower = lower(1:end-1,:);
upper = upper(1:end-1,:);

Y = [lower; upper] ;


%         lower.pop()
%         lower
%         lower.append(p)
%         
%         % Build upper hull
%         upper = [];
%         for p in reversed(points):
%             while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
%                 upper.pop()
%                 upper.append(p)
%             end
%         end
%     end
% end
% 
% % Concatenation of the lower and upper hulls gives the convex hull.
% % Last point of each list is omitted because it is repeated at the beginning of the other list.
% return lower[:-1] + upper[:-1]
