function null_space_vectors = null_tol(A, tol)
    % Computes the null space of a matrix A with a given tolerance
    % Inputs:
    %   A: input matrix
    %   tol: tolerance for determining null space
    % Output:
    %   null_space_vectors: matrix whose columns are the basis vectors of the null space
    
    % Perform SVD on A
    [U, S, V] = svd(A);
    
    % Determine rank of A
    rank_A = sum(diag(S) > tol);
    
    % Extract null space vectors
    null_space_vectors = V(:, rank_A+1:end);
    
    % Remove numerical noise
    null_space_vectors = null_space_vectors(:, abs(diag(null_space_vectors'*null_space_vectors)) > tol);
end