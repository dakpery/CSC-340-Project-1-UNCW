
public class Matrix 
{

	private double[][] matrix;
	private int row;
	private int column; 
	
	

	public Matrix(int r, int c) {
		row = r; column = c;
		matrix = new double[row][column];	
	}
	
	
	public int getRows() {
		return matrix.length;
	}
	
	
	public int getColumns() {

		return matrix[0].length;
	}
	
	
	public Matrix(double[][] matrixArray) {		
		matrix = matrixArray; 
	}
	
	
	public double getValue(int i, int j) {
		return matrix[i][j];
	}
	
	
	public String toString() {
		String printString = (matrix.length + "X" + matrix[0].length + "\n");
		row = matrix.length;
		column = matrix[0].length;
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				printString = printString + matrix[i][j] + "\t\t";
			}	
			printString = printString + "\n";
		}
		return printString;
	}

	
	public Boolean equals(Matrix that) {
		boolean result = false;
		for (int i = 0; i < getRows(); i++) {
			for (int j = 0; j < getColumns(); j++) {
				if (matrix[i][j] == that.matrix[i][j]) {
					result = matrix[i][j] == that.matrix[i][j];
				}
				else {
					return result = false;
				}
			}
		}
		return result;
	}
	
	
	public Matrix transpose() {

		double[][] transArray = new double[getColumns()][getRows()];
		
		for(int i = 0; i < getRows(); i++) {
			for (int j = 0; j < getColumns(); j++) {
				transArray[j][i] = matrix[i][j];
			}
		}
		Matrix newMatrix = new Matrix(transArray);
		return newMatrix; 
	}
	 

	public Matrix add(Matrix that) {

		if (getRows() != that.getRows() || getColumns() != that.getColumns()) {
			System.out.println("error");
			return null; 
		}
		
		double[][] addMatrix = new double[getRows()][getColumns()];
		for (int i = 0; i < getRows(); i++) {
			for (int j = 0; j < getColumns(); j++) {			
				addMatrix[i][j] = matrix[i][j] + that.matrix[i][j];
			}
		}
		return new Matrix(addMatrix);
	}
	
	
	public Matrix subtract (Matrix other) {
		if (getRows() != other.getRows() || getColumns() != other.getColumns()) {
			return null;
		}
		double[][] temp = new double[getRows()][getColumns()];
		for (int i=0; i < getRows(); i++) {
			for (int j=0; j < getColumns(); j++) {
				temp[i][j] = matrix[i][j] - other.matrix[i][j];
			}
		}
		return new Matrix(temp);
	}

	
	public Matrix multiply(double multiplier) {
		Matrix multMatrix = new Matrix(matrix);
		for (int i = 0; i < getRows(); i++) {
			for (int j = 0; j < getColumns(); j++) {
				multMatrix.matrix[i][j] *= multiplier;
			}
		}
		return multMatrix;
	}
	
	
	public Matrix multiply(Matrix other) {
		if (other.getRows() != getColumns())
			return null;

		double[][] c = new double[getRows()][other.getColumns()];
		for (int i=0; i < getRows(); i++) {
			for (int j=0; j < other.getColumns(); j++) {
				for (int k=0; k < getColumns(); k++) {
					c[i][j] = c[i][j] + (matrix[i][k] * other.getValue(k, j));
				}
			}
		}
		Matrix temp = new Matrix(c);
		return temp;
	}
	
	
	public void swapRows(int rowOne, int rowTwo) {
		double[] temp = new double[getColumns()];
		for(int i = 0; i < getColumns(); i++) {
			temp[i] = getValue(rowOne,i);
			matrix[rowOne][i] = matrix[rowTwo][i];
			matrix[rowTwo][i] = temp[i];
		}				
	}
	
	
	public void rowMultiply(int rowNum, double multiplier) {
		for(int i = 0; i < getColumns(); i++) {
			matrix[rowNum][i] *= multiplier; 
		}
	}
	
	
	public void addRows(int rowOne, int rowTwo) {
		for(int i = 0; i < getColumns(); i++) {
			matrix[rowOne][i] += matrix[rowTwo][i];
		}
	}
	
	
	public int findPivot(int j, int pivot, Matrix mat) {
		for(int i = 0; i < mat.getRows(); i++) {
			if (Math.abs(mat.matrix[pivot][j]) < Math.abs(mat.matrix[i][j])) {
				pivot = i;
			} 
		}
		return pivot; 
	}
	
	
	public Matrix augment(Matrix other) {
		double tempMatrix[][] = new double[getRows()][getColumns()+other.getColumns()];
		for(int i = 0; i < tempMatrix.length; i++) {
			for(int j = 0; j < tempMatrix[0].length; j++) {
				if(j < getColumns()) {
					tempMatrix[i][j] = matrix[i][j];
				} else {
					tempMatrix[i][j] = other.matrix[i][j-getColumns()];
				}
			}
		}
		return new Matrix(tempMatrix);
	}
	
	
	public Matrix gaussJordanElimination (Matrix that) {
		int e = 1;
		int pivot;

		// augment current matrix with the values
		Matrix C = augment(that);

		// for every row in the matrix
		// holds the row where we want the diagonal 
		for (int j = 0; j < C.getRows(); j++) {
			
			pivot = j;
			
			findPivot(j,pivot,C);

			if (C.matrix[pivot][j] == 0) {
				e = 0;
			}

			if (pivot > j) {
				C.swapRows(j, pivot);
			}

			C.rowMultiply(j, (1.0/C.getValue(j, j)));

			// for every row 
			// iteratates the rest of the rows 
			for (int i = 0; i < C.getRows(); i++) {
				if (i != j) {
					// for every column
					double Cij = C.matrix[i][j];
					for (int k = 0; k < C.getColumns(); k++) {
						C.matrix[i][k] = C.matrix[i][k] - C.matrix[j][k] * Cij;
					}
				}
			}

		}
		if (e == 0) {
			
			System.out.println("ugh");
			return null;
		}
		else
			return C;
	}
	
	
	public double[] gaussianElim  (Matrix b) {
		int E = 1;
		int p;

		// augment current matrix with the values
		Matrix C = augment(b);

		for (int j=0; j < C.getRows(); j++) {
			p = j;

			// loop through each row and find the pivot
			for (int i=0; i < C.getRows(); i++) {

				// look at the absolute value of the pivot and see if it is larger than current pivot
				// if so, make that the new pivot
				if (Math.abs(C.matrix[p][j]) < Math.abs(C.matrix[i][j])) {
					p = i;
				}
			}

			// if an entire column is zeroes, return 0. no unique soln
			if (C.matrix[p][j] == 0) {
				E = 0;
			}

			// interchange if the pivot is below row j
			if (p > j) {
				C.swapRows(j, p);
			}

			// divide the row by the leading coefficient
			C.rowMultiply(j, (1.0/C.getValue(j, j)));

			for (int i=0; i < C.getRows(); i++) {
				if (i > j) {
					double Cij = C.matrix[i][j];
					double Cjj = C.matrix[j][j];
					for (int k=0; k < C.getColumns(); k++) {
						C.matrix[i][k] = C.matrix[i][k] - (C.matrix[j][k] * (Cij/Cjj));
					}
				}
			}

		}
   
   		// return null if no unique soln
   		if (E == 0)
   			return null;

   		// create partitions D and e
   		// D is the matrix from 0 -> end of original system matrix
		double[][] deaugmentedCoefficients = new double[C.getRows()][C.getColumns() - 1];
		
		// e is the solution matrix 
		double[][] deaugmentedSolutions = new double[C.getRows()][1];
		
		// For j < number of columns in degaugmented matrix D, fill with value c[i][j]
		for (int i = 0; i < C.getRows(); i++) {
			for (int j=0; j < C.getColumns() - 1; j++ ) {
				deaugmentedCoefficients[i][j] = C.getValue(i, j);
			}
		}

		// Since the value of e will only ever be one column
		// we only use one for loop to fill the paritions
		for (int i=0; i < C.getRows(); i++)
			deaugmentedSolutions[i][0] = C.getValue(i, C.getColumns() - 1);
		
		
		
		// solution array
		double[] solution = new double[C.getRows()];
		double sum;
		
		// back substitution
		for (int j=C.getRows() - 1; j >= 0; j--) {
			sum = 0;
			for (int i=j+1; i < C.getRows(); i++) {
				sum = sum + (deaugmentedCoefficients[j][i] * solution[i]);
			}
			solution[j] = (deaugmentedSolutions[j][0] - sum)/deaugmentedCoefficients[j][j];
		}
		
		return solution; 
	}
	
	
	public double findDeterminant() {

		
		if(getColumns() != getRows()) {
			return 0; 
		}
		int pivot; 
		int r = 0; 
		double determinant = 0; 
		
		double[][] currentArray = new double[getRows()][getColumns()];
		for(int i = 0; i < currentArray.length; i++) {
			for(int j = 0; j < currentArray[0].length; j++) {
				currentArray[i][j] = matrix[i][j];
			}
		}
		
		
		Matrix C = new Matrix(currentArray);
		
		
		for(int j = 0; j < C.getRows(); j++) {
			pivot = j; 
			
			findPivot(j,pivot,C);			
			
			if (pivot > j) {
				C.swapRows(j, pivot);
				
				// interchanging rows of the matrix negates the determinant 
				r++;
			}
			
			////////////////////////////////////////////////////
			//
			// Above here is the standard way to find a pivot point 
			//
			////////////////////////////////////////////////////
			//
			// So during backwards substituition
			// 1) For every column in the matrix 
			// 2) Capture each value that is not a zero and not on the diagnol 
			// 3) Apply the forumla: 
			//     a) Cik = Cik - (Cjk * (Cij / Cjj))
			//     b) this formula will allow for the value of Cik to be zeroed if it is in the column J 
			//     c) this also will change the values of any values that are in the row i
			//
			////////////////////////////////////////////////////
			//
			// Below here is iterating the rows with i 
			// Then using another iterator to index into the matrix k 
			//
			////////////////////////////////////////////////////

			
			for (int i=0; i < C.getRows(); i++) {
				if (i > j) {
					
					double Cij = C.matrix[i][j];
					double Cjj = C.matrix[j][j];
					
					for (int k=0; k < C.getColumns(); k++) {				
						C.matrix[i][k] = C.matrix[i][k] - (C.matrix[j][k] * (Cij/Cjj));
					}
				}
			}

		}
		
		
		determinant = C.matrix[0][0];
		
		for (int i=1; i < C.getRows(); i++) {
			determinant = determinant * C.matrix[i][i];

		}
		
		// negates the negation of the determinant if need be. 
		determinant = determinant * Math.pow((-1), r);


		return determinant;	
	}
	

	public Matrix inverse() {
		double[][] tempMatrix = new double[getRows()][getColumns()];
		
		int E = 1;
		int pivot; 
		
		if(getRows() != getColumns()) {
			System.out.println("Not a square matrix");
		}
	
		// Calculating the identity matrix 
		for(int i = 0; i < getColumns(); i++) {
			for(int j = 0; j < getRows(); j++) {
				if(i == j) {
					 tempMatrix[i][j] = 1; 
				}
			}
		}
	
		Matrix identity = new Matrix(tempMatrix);
		
		
		
		Matrix C = augment(identity);
		
		
		for(int j = 0; j < getRows(); j++) {
			pivot = j; 
			
			findPivot(pivot,j,C);
			
			if(C.matrix[pivot][j] == 0) {
				E = 0; 
			}
			
			
			if(pivot > j) {
				C.swapRows(pivot, j);
			}
			
			// divide row j by Cjj 
			C.rowMultiply(j, (1.0/C.getValue(j, j)));

			if(E == 0) {
				return null; 
			}
			

			
			for (int i = 0; i < C.getRows(); i++) {
				if (i != j) {
					double Cij = C.matrix[i][j];					
					for (int k = 0; k < C.getColumns(); k++) {
						C.matrix[i][k] = C.matrix[i][k] - C.matrix[j][k] * Cij;
						
					}
				}
			}
		}
		
		double tempArray [][]  = new double[C.getRows()][C.getColumns() / 2];
			
		
		for(int i = 0; i < C.getRows(); i++) {
			for(int j = C.getColumns() / 2;  j < C.getColumns(); j++) {
				tempArray[i][j - C.getColumns() / 2] = C.matrix[i][j];
			}
		}
		
		
		
		return new Matrix(tempArray);		
	}
	
	public double trace() {
		if (getColumns() != getRows())
			return 0.0; 
		double sum = 0.0; 
		for(int i = 0; i < getRows(); i++)
			sum += matrix[i][i];
		return sum; 
	}
}







