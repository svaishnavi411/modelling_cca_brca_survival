function [U, V] = sparse_cca_pd(X, Y, paras)

[U, V] = SparseCCAGroupMultiple(X, Y, paras.lambda1, paras.lambda2, ...
                                   [], [], [], [], paras.dim, 'pd');