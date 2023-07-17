function [U, V] = sparse_cca_hd(X, Y, paras)

[U, V] = SparseCCAGroupMultiple(X, Y, paras.lambda1, paras.lambda2, ...
                                   [], [], [], [], paras.dim, 'hd');