pts=[[90  34  89];[21   3  78];[76  35  74];[28   4   7];[65  60  22];[59  92   5];[84   0  32]];

tri = DelaunayTri(pts);
tetrahedra = tri.Triangulation;
centers = tri.circumcenters;

num_tetrahedra = size(tetrahedra, 1)
for t=1:num_tetrahedra
    tetrahedron = tetrahedra(t,:);
    for i=1:4
        pt_index = tetrahedron(i);
        pt = pts(pt_index,:);
        fprintf('(');
        fprintf('%.0f ', pt);
        fprintf('\b) ');
    end
    fprintf(' - (');
    fprintf('%f ', centers(t,:));
    fprintf('\b)\n');
end
