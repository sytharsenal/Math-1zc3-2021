
Eigenvalue & eigenvectors
    e = eig(A)    #eigenvalue
    [V,D] = eig(A)    #eigenvectors

Transpose
    transpose(A) or [A]'

Cross-product of u,v
    cross(u,v)

Generating sequence
    x:y #from x to y
    x:y:z #starting as x increase by y ending at z

Cross-product
    C = cross(A,B)

Repeat copies of an array
    repmat(A,n)

Lower triangular part of a matrix
    tril(A)
    #returns the lower triangular portion of matrix A.
    tril(A,k)
    #returns the elements on and below the kth diagonal of A.

Gram-Schmidt
    for j =1:n 
        v = A(:,j);
        # v begins as column j of A
        for i = 1:j-1
            R(i,j)= Q(:i)'*A(:j);
            # modify A(;,j) to v for more accuracy
            v=v-R(i,j)*Q(:,i);
            #subtract the projection 
        end
        R(j,j)=norm(v);
        #v perpendicular to all of q_1, ..., q_j-1
        Q(:,j)=v/R(j,j); 
        #normalize v to be the next unit vector q_j
    end

Calculator orthonormal basis for range of A
    A = [1 0 1;-1 -2 0; 0 1 -1];
    r = rank(A)

    Q = orth(A)

Projecting A to B
    A=[5,-1,1];
    B=[2,-1,-2];
    #calculation of the projection of A into B
    C=(sum(A.*B)/(norm(B)^2))*B;

Column space of A
    colspace(A)

Row space of A
    A_tr = transpose(A) 
    A_rs = colspace(A_tr)
    #or ------------------
    colspace(A')

Null space of A
    null(A,'r')
    null(sym(A))
    null(A)  %reduced decimal form