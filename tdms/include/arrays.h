/**
 * @file arrays.h
 * @brief Classes describing arrays, vertices etc.
 */
#pragma once

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include <fftw3.h>

#include "globals.h"
#include "matlabio.h"
#include "utils.h"

template<typename T>
class XYZTensor3D {
public:
  T ***x = nullptr;
  T ***y = nullptr;
  T ***z = nullptr;

  T ***operator[](char c) const {
    switch (c) {
      case 'x':
        return x;
      case 'y':
        return y;
      case 'z':
        return z;
      default:
        throw std::runtime_error("Have no element " + std::string(1, c));
    }
  }
  T ***operator[](AxialDirection d) const {
    switch (d) {
      case AxialDirection::X:
        return x;
      case AxialDirection::Y:
        return y;
      case AxialDirection::Z:
        return z;
      default:
        throw std::runtime_error("Have no element " + std::to_string(d));
    }
  }

  /**
   * @brief Allocates x, y, and z as (K_total+1) * (J_total+1) * (I_total+1)
   * arrays
   *
   * @param I_total,J_total,K_total Dimensions of the tensor size to set
   */
  void allocate(int I_total, int J_total, int K_total) {
    x = (T ***) malloc(K_total * sizeof(T **));
    y = (T ***) malloc(K_total * sizeof(T **));
    z = (T ***) malloc(K_total * sizeof(T **));
    for (int k = 0; k < K_total; k++) {
      x[k] = (T **) malloc(J_total * sizeof(T *));
      y[k] = (T **) malloc(J_total * sizeof(T *));
      z[k] = (T **) malloc(J_total * sizeof(T *));
      for (int j = 0; j < J_total; j++) {
        x[k][j] = (T *) malloc(I_total * sizeof(T));
        y[k][j] = (T *) malloc(I_total * sizeof(T));
        z[k][j] = (T *) malloc(I_total * sizeof(T));
      }
    }
  }
};

class XYZVectors {
public:
  double *x = nullptr;
  double *y = nullptr;
  double *z = nullptr;

  /**
   * Default constructor
   */
  XYZVectors() = default;

  /**
   * Set the pointer for one of the vectors in this collection with a name of c
   * @param c Character labeling the vector
   * @param ptr Pointer to assign
   */
  void set_ptr(char c, double *ptr);
  /**
   * Set the pointer for one of the vectors in this collection with a name of c
   * @param d AxialDirection labeling the vector
   * @param ptr Pointer to assign
   */
  void set_ptr(AxialDirection d, double *ptr);

  /**
   * @brief Determines whether all elements in the x, y, or z vector are less
   * than a given value.
   *
   * @param comparison_value Value to compare elements to
   * @param vector_length Number of elements to compare against
   * @param component Vector to compare elements against; x, y, or z
   * @param buffer_start Only compare elements between buffer_start (inclusive)
   * and buffer_start+vector_length-1 (inclusive)
   * @return true All elements are less than the comparison_value
   * @return false At least one element is not less than the comparison_value
   */
  bool all_elements_less_than(double comparison_value, int vector_length,
                              AxialDirection component,
                              int buffer_start = 0) const;
  /**
   * @brief Determines whether all elements in the x, y, AND z vectors are less
   * than a given value.
   *
   * @param comparison_value Value to compare elements to
   * @param nx,ny,nz Number of elements in the nx, ny, and nz vectors
   * respectively
   * @return true All elements are less than the comparison_value
   * @return false At least one element is not less than the comparison_value
   */
  bool all_elements_less_than(double comparison_value, int nx, int ny,
                              int nz) const;
};

// TODO: docstring
class MaterialCollection {
protected:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                               const std::string &prefix);
};

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwells equations for E fields. The symbol chosen in
 * the original reference is \f$C\f$.
 *
 * @details Algebraic terms \f$C_{a,b,c}\f$ defined in Section 4.2 of Munro, P,.
 * "Application of numerical methods to high numerical aperture imaging", 2006,
 * PhD thesis, Imperial College London.
 *
 * The definitions are equations 4.13, 4.14 (pp 82-3). Part of Maxwell's E-field
 * equations in equations 4.7-9.
 */
class CCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
  XYZVectors c;
};

/*! @copydoc CCollectionBase */
class CCollection : public CCollectionBase {
private:
  void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                        const std::string &prefix);

public:
  bool is_multilayer = false;
  bool is_disp_ml = false;

  explicit CCollection(const mxArray *ptr);
};

/*! @copydoc CCollectionBase */
class CMaterial : public CCollectionBase, MaterialCollection {
public:
  explicit CMaterial(const mxArray *ptr);
};

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwells equations for H fields. The symbol chosen in
 * the original reference is \f$D\f$.
 *
 * @details Algebraic terms \f$D_{a,b}\f$ defined in Section 4.2 of Munro, P,.
 * "Application of numerical methods to high numerical aperture imaging", 2006,
 * PhD thesis, Imperial College London.
 *
 * The definitions are equations 4.15, 4.16 (pp 82-3). Part of Maxwell's H-field
 * equations in equations 4.10-12.
 */
class DCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
};

/*! @copydoc DCollectionBase */
class DCollection : public DCollectionBase {
private:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                               const std::string &prefix);

public:
  explicit DCollection(const mxArray *ptr);
};

/*! @copydoc DCollectionBase */
class DMaterial : public DCollectionBase, MaterialCollection {
public:
  explicit DMaterial(const mxArray *ptr);
};

class DispersiveMultiLayer {
public:
  double *alpha = nullptr;
  double *beta = nullptr;
  double *gamma = nullptr;
  XYZVectors kappa;
  XYZVectors sigma;

public:
  explicit DispersiveMultiLayer(const mxArray *ptr);

  /**
   * @brief Determines whether the (background) medium is dispersive
   *
   * @param K_tot Number of Yee cells in the z-direction (number of entries in
   * this->gamma)
   * @param near_zero_tolerance Tolerance for non-zero gamma (attenuation)
   * values
   * @return true Background is dispersive
   * @return false Background is not dispersive
   */
  bool is_dispersive(int K_tot, double near_zero_tolerance = 1e-15);
};

template<typename T>
class Matrix {
protected:
  int n_rows = 0;
  int n_cols = 0;
  T **matrix = nullptr;

public:
  /**
   * @brief Construct a new Matrix object, without assigned elements
   *
   */
  Matrix() = default;
  /**
   * @brief Construct a new Matrix object, providing the dimensions
   *
   * @param n_rows,n_cols Number of rows and columns in the matrix
   * @param initial_value The initial value of the elements, defaults to 0 to
   * avoid initalised but unassigned values
   */
  Matrix(int n_rows, int n_cols) { allocate(n_rows, n_cols); }

  inline T *operator[](int value) const { return matrix[value]; }
  /**
   * @brief Check whether this matrix has elements assigned
   *
   * @return true If this matrix has assigned elements
   * @return false This matrix is currently unassigned
   */
  bool has_elements() { return matrix != nullptr; };

  /**
   * Allocate the memory for this matrix
   *
   * @param n_rows Number of rows
   * @param n_cols Number of columns
   */
  void allocate(int n_rows, int n_cols, T initial_value = 0) {
    this->n_rows = n_rows;
    this->n_cols = n_cols;

    matrix = (T **) malloc(sizeof(T *) * n_rows);
    for (int i = 0; i < n_rows; i++) {
      matrix[i] = (T *) malloc(sizeof(T) * n_cols);
    }
  };

  int get_n_cols() const { return n_cols; }
  int get_n_rows() const { return n_rows; }

  /**
   * Destructor. Must be defined in the header
   */
  ~Matrix() {
    if (has_elements()) {
      for (int i = 0; i < n_rows; i++) { free(matrix[i]); }
      free(matrix);
    }
  };
};

class GratingStructure : public Matrix<int> {

public:
  GratingStructure(const mxArray *ptr, int I_tot);

  ~GratingStructure();
};

template<typename T>
class Vector {
protected:
  int n = 0;          // Number of elements
  T *vector = nullptr;// Internal array

public:
  Vector() = default;

  explicit Vector(const mxArray *ptr) {
    n = (int) mxGetNumberOfElements(ptr);
    vector = (T *) malloc((unsigned) (n * sizeof(T)));

    auto matlab_ptr = mxGetPr(ptr);
    for (int i = 0; i < n; i++) { vector[i] = (T) matlab_ptr[i]; }
  }

  bool has_elements() { return vector != nullptr; }

  inline T operator[](int value) const { return vector[value]; };

  inline int size() const { return n; };
};

class FrequencyExtractVector : public Vector<double> {
public:
  FrequencyExtractVector(const mxArray *ptr, double omega_an);

  double max();
};

struct FrequencyVectors {
  std::vector<double> x;
  std::vector<double> y;
};

/**
 * @brief Defines the numerical aperture of the objective, assuming that the
 * lens is centred on the origin of the PSTD simulation.
 *
 * In particular, since the fibre modes are imaged onto a Fourier plane of both
 * the physical fibre and the sample, the field scattered by the sample and
 * collected by the objective lens can have only a finite spatial support in the
 * aperature of the objective lens.
 *
 * Pupil[j][i] thus takes the value 1 for those (i,j) indices (note the order
 * swapping) within the aperture of the lens.
 */
class Pupil : public Matrix<double> {
public:
  Pupil() = default;

  void initialise(const mxArray *ptr, int n_rows, int n_cols);

  ~Pupil();
};

template<typename T>
class Tensor3D {
protected:
  int n_layers = 0;
  int n_cols = 0;
  int n_rows = 0;
  T ***tensor = nullptr;

public:
  bool is_matlab_initialised = false;

  Tensor3D() = default;

  Tensor3D(T ***tensor, int n_layers, int n_cols, int n_rows) {
    initialise(tensor, n_layers, n_cols, n_rows);
  }

  void initialise(T ***_tensor, int _n_layers, int _n_cols, int _n_rows) {
    tensor = _tensor;
    n_layers = _n_layers;
    n_cols = _n_cols;
    n_rows = _n_rows;
  }

  inline T **operator[](int value) const { return tensor[value]; };

  bool has_elements() { return tensor != nullptr; };

  void zero() {
    for (int k = 0; k < n_layers; k++)
      for (int j = 0; j < n_cols; j++)
        for (int i = 0; i < n_rows; i++) { tensor[k][j][i] = 0; }
  }

  void allocate(int nK, int nJ, int nI) {
    n_layers = nK, n_cols = nJ, n_rows = nI;
    tensor = (T ***) malloc(n_layers * sizeof(T **));

    for (int k = 0; k < n_layers; k++) {
      tensor[k] = (T **) malloc(n_cols * sizeof(T *));
    }

    for (int k = 0; k < n_layers; k++) {
      for (int j = 0; j < n_cols; j++) {
        tensor[k][j] = (T *) malloc(n_rows * sizeof(T));
      }
    }
  };

  /**
   * @brief Computes the Frobenius norm of the tensor
   *
   * fro_norm = \f$\sqrt{ \sum_{i=0}^{I_tot}\sum_{j=0}^{J_tot}\sum_{k=0}^{K_tot}
   * |t[k][j][i]|^2 }\f$
   */
  double frobenius() {
    T norm_val = 0;
    for (int i1 = 0; i1 < n_layers; i1++) {
      for (int i2 = 0; i2 < n_cols; i2++) {
        for (int i3 = 0; i3 < n_rows; i3++) {
          norm_val += abs(tensor[i1][i2][i3]) * abs(tensor[i1][i2][i3]);
        }
      }
    }
    return sqrt(norm_val);
  }

  ~Tensor3D() {
    if (tensor == nullptr) return;
    if (is_matlab_initialised) {
      free_cast_matlab_3D_array(tensor, n_layers);
    } else {
      for (int k = 0; k < n_layers; k++) {
        for (int j = 0; j < n_cols; j++) { free(tensor[k][j]); }
        free(tensor[k]);
      }
      free(tensor);
    }
  }
};

/**
 * @brief Stores the fibre modes in the Fourier plane of the objective lens.
 *
 * The "Tilde" indicates that these quantities are in a Fourier plane relative
 * to where the optical fibre is actually located, meaning that is has a Fourier
 * relationship relative to the physical fibre mode(s).
 */
class DTilde {
protected:
  int n_det_modes = 0;//< Number of modes specified
  static void set_component(Tensor3D<std::complex<double>> &tensor,
                            const mxArray *ptr, const std::string &name,
                            int n_rows, int n_cols);

public:
  /** @brief Fetch the number of modes */
  inline int num_det_modes() const { return n_det_modes; };

  /*! 3-dimensional vector of fibre modes indexed by (j, i, i_m).
   * i and j index over the x and y plane respectively.
   * i_m indexes the different modes specified in the input file.*/
  Tensor3D<std::complex<double>> x;
  /*! @copydoc x */
  Tensor3D<std::complex<double>> y;

  void initialise(const mxArray *ptr, int n_rows, int n_cols);
};

class IncidentField {
protected:
  void set_component(Tensor3D<double> &component, const mxArray *ptr,
                     const std::string &name);

public:
  Tensor3D<double> x;
  Tensor3D<double> y;

  explicit IncidentField(const mxArray *ptr);
};

/**
 * List of field components as integers
 */
class FieldComponentsVector : public Vector<int> {
public:
  FieldComponentsVector() = default;

  void initialise(const mxArray *ptr);

  /**
   * Get the index of a particular integer in this vector. If it does not exist
   * then return -1. Returns the first occurrence.
   * @param value value to find in the vector
   * @return index or -1
   */
  int index(int value);
};

class Vertices : public Matrix<int> {
public:
  Vertices() = default;

  void initialise(const mxArray *ptr);

  int n_vertices() { return n_rows; }

  ~Vertices() {
    if (has_elements()) { free_cast_matlab_2D_array(matrix); }
    matrix = nullptr;
  };
};

class DetectorSensitivityArrays {
public:
  fftw_complex *v = nullptr;          // Flat fftw vector
  fftw_plan plan = nullptr;           // fftw plan for the setup
  std::complex<double> **cm = nullptr;// Column major matrix

  void initialise(int n_rows, int n_cols);

  ~DetectorSensitivityArrays();
};

/**
 * Matrix of c coefficients. See the pdf documentation for their definition
 */
class CCoefficientMatrix : public Matrix<double> {};

/**
 * Temporary storage 'vector'
 */
class EHVec : public Matrix<fftw_complex> {
public:
  ~EHVec();
};

/**
 * Container for storing snapshots of the full-field
 */
class FullFieldSnapshot {
public:
  std::complex<double> Ex = 0.;//< x-component of the electric field
  std::complex<double> Ey = 0.;//< y-component of the electric field
  std::complex<double> Ez = 0.;//< z-component of the electric field
  std::complex<double> Hx = 0.;//< x-component of the magnetic field
  std::complex<double> Hy = 0.;//< y-component of the magnetic field
  std::complex<double> Hz = 0.;//< z-component of the magnetic field

  FullFieldSnapshot() = default;

  /**
   * @brief Return the component of the field corresponding to the index
   * provided.
   *
   * 0 = Ex, 1 = Ey, 2 = Ez, 3 = Hx, 4 = Hy, 5 = Hz.
   * This is the indexing order that other storage containers use.
   *
   * Throws an error if provided an index <0 or >5.
   *
   * @param index Field component index to fetch
   * @return std::complex<double> The field component
   */
  std::complex<double> operator[](int index) {
    switch (index) {
      case 0:
        return Ex;
        break;
      case 1:
        return Ey;
        break;
      case 2:
        return Ez;
        break;
      case 3:
        return Hx;
        break;
      case 4:
        return Hy;
        break;
      case 5:
        return Hz;
        break;
      default:
        throw std::runtime_error("Index " + std::to_string(index) +
                                 " does not correspond to a field component.");
        break;
    }
  }

  /**
   * @brief Multiplies the electric field components by `factor`.
   * @param factor to scale the field by
   */
  void multiply_E_by(std::complex<double> factor) {
    Ex *= factor;
    Ey *= factor;
    Ez *= factor;
  }
  /**
   * @brief Multiplies the magnetic field components by `factor`.
   * @param factor to scale the field by
   */
  void multiply_H_by(std::complex<double> factor) {
    Hx *= factor;
    Hy *= factor;
    Hz *= factor;
  }
};
