function y = crosspod(o, a, b)% 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
% Returns a positive value, if OAB makes a counter-clockwise turn,
% negative for clockwise turn, and zero if the points are collinear.

y = (a(1) - o(1)) * (b(2) - o(2)) - (a(2) - o(2)) * (b(1) - o(1));