"""
Derived module from dmdbase.py for classic dmd.
"""
from .dmdbase import DMDBase
import scipy as sp

class DMD(DMDBase):
    """
    Dynamic Mode Decomposition

    :param svd_rank: the rank for the truncation; If 0, the method computes the
        optimal rank and uses it for truncation; if positive interger, the
        method uses the argument for the truncation; if float between 0 and 1,
        the rank is the number of the biggest singular values that are needed
        to reach the 'energy' specified by `svd_rank`; if -1, the method does
        not compute truncation.
    :type svd_rank: int or float
    :param int tlsq_rank: rank truncation computing Total Least Square. Default
        is 0, that means TLSQ is not applied.
    :param bool exact: flag to compute either exact DMD or projected DMD.
        Default is False.
    :param bool opt: flag to compute optimal amplitudes. See :class:`DMDBase`.
        Default is False.
    """
    def __init__(self, alpha, svd_rank=0, tlsq_rank=0, exact=False, opt=False):
        self.svd_rank = svd_rank
        self.tlsq_rank = tlsq_rank
        self.exact = exact
        self.opt = opt
        self.original_time = None
        self.dmd_time = None
		self.alpha = alpha
        self._eigs = None
        self._Atilde = None
        self._modes = None  # Phi
        self._b = None  # amplitudes
        self._snapshots = None
        self._snapshots_shape = None

    @staticmethod
    def _eig_from_lowrank_op(Atilde, Y, U, s, V, exact):
        """
        Private method that computes eigenvalues and eigenvectors of the
        high-dimensional operator from the low-dimensional operator and the
        input matrix.

        :param numpy.ndarray Atilde: the lowrank operator.
        :param numpy.ndarray Y: input matrix Y.
        :param numpy.ndarray U: 2D matrix that contains the left-singular
            vectors of X, stored by column.
        :param numpy.ndarray s: 1D array that contains the singular values of X.
        :param numpy.ndarray V: 2D matrix that contains the right-singular
            vectors of X, stored by row.
        :param bool exact: if True, the exact modes are computed; otherwise,
            the projected ones are computed.
        :return: eigenvalues, eigenvectors
        :rtype: numpy.ndarray, numpy.ndarray
        """
        lowrank_eigenvalues, lowrank_eigenvectors_l, lowrank_eigenvectors_r = sp.linalg.eig(Atilde, left=True, right=True)

        # Compute the eigenvectors of the high-dimensional operator
        if exact:
            eigenvectors = (
                (Y.dot(V) * np.reciprocal(s)).dot(lowrank_eigenvectors_l))
			raise('currently adjoint DMD use only projected modes')

        else:
            eigenvectors = U.dot(lowrank_eigenvectors_l)
			eigenvectors_adj = U.dot(lowrank_eigenvectors_r.conj().T)
		
        # The eigenvalues are the same
        eigenvalues = lowrank_eigenvalues

        return eigenvalues, eigenvectors,eigenvectors_adj
    

	@staticmethod  
	def mode_projection_coefficients(self):
		pass
	@staticmethod  
	def projected_dynamics(self,X,t, alpha, L):
		pass







	def fit(self, X):
        """
        Compute the Dynamic Modes Decomposition to the input data.

        :param X: the input snapshots.
        :type X: numpy.ndarray or iterable
        """
        self._snapshots, self._snapshots_shape = self._col_major_2darray(X)

        n_samples = self._snapshots.shape[1]
        X = self._snapshots[:, :-1]
        Y = self._snapshots[:, 1:]

        X, Y = self._compute_tlsq(X, Y, self.tlsq_rank)

        U, s, V = self._compute_svd(X, self.svd_rank)

        self._Atilde = self._build_lowrank_op(U, s, V, Y)

        self._eigs, self._modes, self._adjmodes = self._eig_from_lowrank_op(
            self._Atilde, Y, U, s, V, self.exact)
		
		L = mode_projection_coefficients() # r by r
		a = projected_dynamics(self,X,t, alpha, L) # a is a matrix r by t 
        
		# Default timesteps
        self.original_time = {'t0': 0, 'tend': n_samples - 1, 'dt': 1}
        self.dmd_time = {'t0': 0, 'tend': n_samples - 1, 'dt': 1}
        
        self._b = self._compute_amplitudes(self._modes, self._snapshots,
                                           self._eigs, self.opt)

        return self


    @property
    def reconstructed_data(self):
		L = mode_projection_coefficients() # r by r
		a = projected_dynamics(self,self.X, t, alpha, L) # a is a matrix r by t 
		return self.modes @ a




