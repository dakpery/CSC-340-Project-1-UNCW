import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

public class Test {
	
	private static ArrayList<Matrix> class1 = new ArrayList<Matrix>();
	private static ArrayList<Matrix> class2 = new ArrayList<Matrix>(); 
	
	private static ArrayList<Matrix> class1Boundaries = new ArrayList<Matrix>();
	private static ArrayList<Matrix> class2Boundaries = new ArrayList<Matrix>();

	
	public static void main(String args[]) throws FileNotFoundException {

		Test t1 = new Test();
		t1.readFile();
		
		// Mean Vectors
		Matrix meanVector1 = t1.findClass1MeanVector(class1);
		Matrix meanVector2 = t1.findClass2MeanVector(class2);
		System.out.println("meanVector 1 : " + meanVector1.toString());
		System.out.println("meanVector 2 : " + meanVector2.toString());
		
		// Covariance Matrices
		Matrix covariance1 = t1.findCovarianceMatrix1(class1, meanVector1);
		Matrix covariance2 = t1.findCovarianceMatrix2(class2, meanVector2);
		System.out.println("Covariance 1 Matrix: " + covariance1.toString());
		System.out.println("Covariance 1 Matrix: " + covariance2.toString());

		// Determinants of the covariance matrices		
		double determinant1 = covariance1.findDeterminant(); 
		double determinant2 = covariance2.findDeterminant(); 
		System.out.println("Covariance 1 Determinant: " + determinant1);
		System.out.println("Covariance 2 Determinant: " + determinant2);
//		
//		System.out.println("");
//		
		//Inverses for the covariance matrices 
		Matrix covInverse1 = covariance1.inverse(); 
		Matrix covInverse2 = covariance2.inverse(); 
		System.out.println("Covariance 1 Inverse: " + covInverse1.toString());
		System.out.println("Covariance 1 inverse: " + covInverse2.toString());
//
		System.out.println(t1.classifier1(meanVector1, covInverse1, meanVector2, covInverse2, determinant1));
		System.out.println(t1.classifier2(meanVector1, covInverse1, meanVector2, covInverse2, determinant2));
		
		System.out.println("");
		System.out.println("##############################################################");
		System.out.println("");
		
		ArrayList<Matrix> misClassifieds1 = t1.findMisClassifieds1(meanVector1, meanVector2, determinant1, determinant2, covInverse1, covInverse2);
		for(int i = 0; i < misClassifieds1.size(); i++) {
			System.out.println("Class 1 misclassifieds\n" + "\t" + misClassifieds1.get(i).toString());
		}
		
		System.out.println("");
		System.out.println("##############################################################");
		System.out.println("");

		ArrayList<Matrix> misClassifieds2 = t1.findMisClassifieds2(meanVector1, meanVector2, determinant1, determinant2, covInverse1, covInverse2);
		for(int i = 0; i < misClassifieds2.size(); i++) {
			System.out.println("Class 2 misclassified\n" + "\t" + misClassifieds2.get(i).toString());
		}
		
		System.out.println("");
		System.out.println("##############################################################");
		System.out.println("");
		
		for(int i = 0; i < class1Boundaries.size(); i++) {
			System.out.println("Class 1 Boundaries\n" + class1Boundaries.get(i).toString());
		}
		
		System.out.println("");
		System.out.println("##############################################################");
		System.out.println("");
		
		for(int i = 0; i < class2Boundaries.size(); i++) {
			System.out.println("Class 2 Boundaries\n" + class2Boundaries.get(i).toString());
		}
		
		System.out.println("");
		System.out.println("##############################################################");
		System.out.println("");
		

		
		Matrix systemOfEquations = t1.systemOfEquations(); // augmented gauss jordan system 
		System.out.println("Gauss Jordan \n" + systemOfEquations.toString());
		
		t1.printNumberNine();
		t1.printConditionNumber(); 
		
		
		

		


	}
	
	public void readFile() throws FileNotFoundException {
		FileReader inputFile = new FileReader("Resources/data.txt");
		Scanner fileReader = new Scanner(inputFile);
		
		while(fileReader.hasNext()) {
			String line = fileReader.nextLine(); 
			Scanner stringReader = new Scanner(line);
			double classMatrix1[][] = new double[2][1];
			double classMatrix2[][] = new double[2][1];
			classMatrix1[0][0] = stringReader.nextDouble();
			classMatrix1[1][0] = stringReader.nextDouble(); 
			class1.add(new Matrix(classMatrix1));
			classMatrix2[0][0] = stringReader.nextDouble(); 
			classMatrix2[1][0] = stringReader.nextDouble();
			class2.add(new Matrix(classMatrix2));
		}
		fileReader.close();
	}
	
	public Matrix findClass1MeanVector(ArrayList<Matrix> class1) {
		Matrix class1MeanVector = class1.get(0);
		for(int i = 1; i < class1.size(); i++) {
			class1MeanVector = class1MeanVector.add(class1.get(i));
		}
		double class1Scaler =  1.0 / class1.size();
		class1MeanVector = class1MeanVector.multiply(class1Scaler);

		return class1MeanVector; 
	}
	
	public Matrix findClass2MeanVector(ArrayList<Matrix> class1) {
		Matrix class2MeanVector = class1.get(0);
		for(int i = 1; i < class1.size(); i++) {
			class2MeanVector = class2MeanVector.add(class1.get(i));
		}
		double class1Scaler =  1.0 / class1.size();
		class2MeanVector = class2MeanVector.multiply(class1Scaler);

		return class2MeanVector; 
	}
	
	public Matrix findCovarianceMatrix1(ArrayList<Matrix> class1, Matrix meanVector) {
		Matrix covariance1 = class1.get(0).subtract(meanVector);
		covariance1 = covariance1.multiply(covariance1.transpose());
		
		for(int i=1; i < class1.size(); i++) {
			Matrix temp = class1.get(i);
			temp = temp.subtract(meanVector);
			temp = temp.multiply(temp.transpose());
			covariance1 = covariance1.add(temp);
		}

		covariance1 = covariance1.multiply(1.0 / class1.size()); 
		
		return covariance1; 
	}
	
	public Matrix findCovarianceMatrix2(ArrayList<Matrix> class2, Matrix meanVector) {
		Matrix covariance2 = class2.get(0).subtract(meanVector);
		covariance2 = covariance2.multiply(covariance2.transpose());
		
		for(int i=1; i < class2.size(); i++) {
			Matrix temp = class2.get(i);
			temp = temp.subtract(meanVector);
			temp = temp.multiply(temp.transpose());
			covariance2 = covariance2.add(temp);
		}

		covariance2 = covariance2.multiply(1.0 / class1.size()); 
		
		return covariance2;  
	}
	
	public String classifier1(Matrix m1, Matrix cov1, Matrix m2, Matrix cov2, double det) { 
		Matrix g1 = m1.subtract(m1);
		g1 = g1.transpose(); 
		g1 = g1.multiply(cov1);
		g1 = g1.multiply(m1.subtract(m1));
		g1 = g1.multiply(-.5);
		
		double secondHalf = 0.0; 
		secondHalf = Math.log(Math.abs(det));
		secondHalf += Math.log(.5);
		double g1Value = g1.getValue(0, 0) - secondHalf;
		
		
		Matrix g2 = m1.subtract(m2);
		g2 = g2.transpose(); 
		g2 = g2.multiply(cov2);
		g2 = g2.multiply(m1.subtract(m2));
		g2 = g2.multiply(-.5);
		
		secondHalf = 0.0; 
		secondHalf = Math.log(Math.abs(det));
		secondHalf += Math.log(.5);
		double g2Value = g2.getValue(0, 0) - secondHalf;

		if(g1.getValue(0,0) > g2.getValue(0, 0)) {
			System.out.println("m1 was classified in class 1");
		}
		else {
			System.out.println("m1 was classified in class 2");
		}
		
		return ("Class 1 = " + g1Value + " " + g2Value);
	}
	
	public String classifier2(Matrix m1, Matrix cov1, Matrix m2, Matrix cov2, double det) { 
		Matrix g1 = m2.subtract(m2);
		g1 = g1.transpose(); 
		g1 = g1.multiply(cov2);
		g1 = g1.multiply(m2.subtract(m2));
		g1 = g1.multiply(-.5);
		
		double secondHalf = 0.0; 
		secondHalf = Math.log(Math.abs(det));
		secondHalf += Math.log(.5);
		double g1Value = g1.getValue(0, 0) - secondHalf;
		
		
		Matrix g2 = m2.subtract(m1);
		g2 = g2.transpose(); 
		g2 = g2.multiply(cov1);
		g2 = g2.multiply(m2.subtract(m1));
		g2 = g2.multiply(-.5);
		
		secondHalf = 0.0; 
		secondHalf = Math.log(Math.abs(det));
		secondHalf += Math.log(.5);
		double g2Value = g2.getValue(0, 0) - secondHalf;

		if(g1.getValue(0,0) > g2.getValue(0, 0)) {
			System.out.println("\nm2 was classified in class 2");
		}
		else {
			System.out.println("\nm2 was classified in class 1");
		}
		
		return ("Class 1 = " + g1Value + " " + g2Value);
	}
	
	public ArrayList<Matrix> findMisClassifieds1(Matrix m1, Matrix m2,double d1,double d2 ,Matrix cov1, Matrix cov2) {
		System.out.println("");
        Double epsilon = 0.01;
        Double scale = 0.005;
        Double comp = 0.0;
        
        ArrayList<Matrix> class1Outliers = new ArrayList<Matrix>();

		
		Matrix g1; 	
		Matrix g2; 
		double g1Answer;
		double g2Answer;
		double secondHalf; 
		
		for(int i = 0; i < class1.size(); i++) {
			
			g1 = class1.get(i).subtract(m1);
			g1 = g1.multiply(-0.5);
			g1 = g1.transpose();
			g1 = g1.multiply(cov1); // creates a 1 x 2
			g1 = g1.multiply(class1.get(i).subtract(m1));
			g1Answer = g1.getValue(0, 0);
			
			secondHalf = (0.5) * Math.log(d1);
			secondHalf += Math.log(0.5);
			g1Answer -= secondHalf; 
			
			////////////////////////////////////////
			
			secondHalf = 0.0; 
			
			g2 = class1.get(i).subtract(m2);
			g2 = g2.multiply(-0.5);
			g2 = g2.transpose();
			g2 = g2.multiply(cov2); // creates a 1 x 2
			g2 = g2.multiply(class1.get(i).subtract(m2));
			g2Answer = g2.getValue(0, 0);
			
			secondHalf = (0.5) * Math.log(d2);
			secondHalf += Math.log(0.5);
			g2Answer -= secondHalf; 

			
            if (g1Answer < g2Answer) {

            	class1Outliers.add(class1.get(i));

            }
            
            g1Answer *= scale;
            g2Answer *= scale;
            comp = Math.abs(g1Answer - g2Answer);
            if (comp < epsilon)
            {


                this.class1Boundaries.add(class1.get(i));
            }
		}
		

		
		return class1Outliers; 

	}

	public ArrayList<Matrix> findMisClassifieds2(Matrix m1, Matrix m2,double d1,double d2 ,Matrix cov1, Matrix cov2) {
		System.out.println("");
        Double epsilon = 0.01;
        Double scale = 0.005;
        Double comp = 0.0;
        
        ArrayList<Matrix> class2Outliers = new ArrayList<Matrix>();
        ArrayList<Matrix> class2Boundaries = new ArrayList<Matrix>();

		
		Matrix g1; 	
		Matrix g2; 
		double g1Answer;
		double g2Answer;
		double secondHalf; 
		
		for(int i = 0; i < class2.size(); i++) {
			
			g1 = class2.get(i).subtract(m2);
			g1 = g1.multiply(-0.5);
			g1 = g1.transpose();
			g1 = g1.multiply(cov2); // creates a 1 x 2
			g1 = g1.multiply(class2.get(i).subtract(m2));
			g1Answer = g1.getValue(0, 0);
			
			secondHalf = (0.5) * Math.log(d2);
			secondHalf += Math.log(0.5);
			g1Answer -= secondHalf; 
			
			////////////////////////////////////////
			
			secondHalf = 0.0; 
			
			g2 = class2.get(i).subtract(m1);
			g2 = g2.multiply(-0.5);
			g2 = g2.transpose();
			g2 = g2.multiply(cov1); // creates a 1 x 2
			g2 = g2.multiply(class2.get(i).subtract(m1));
			g2Answer = g2.getValue(0, 0);
			
			secondHalf = (0.5) * Math.log(d1);
			secondHalf += Math.log(0.5);
			g2Answer -= secondHalf; 

            if (g1Answer < g2Answer) {

            	class2Outliers.add(class2.get(i));

            }

            
            g1Answer *= scale;
            g2Answer *= scale;
            comp = Math.abs(g1Answer - g2Answer);
            



            
            if (comp < epsilon)
            {
                this.class2Boundaries.add(class2.get(i));
            }
		}
		
		
		return class2Outliers; 

	}
	
	public Matrix systemOfEquations() {
		/**
		 *  2x + y - z - w + a + 0 - c - d = 1
		 *  
		 *  x + 0 +2z + 0 - a -2b + 2c + 2d = -1
		 *  
		 *  0 - 2y + 5z + 4w - a + 0 + 3c + d = 2
		 *  
		 *  x + y -7z + 3w + 2a + b - c + 0 = -2
		 *  
		 *  x + y + 2z + 3w - 2a + 2b + 2c + 9d = 3
		 *  
		 *  0 -3y - 2z + 2w + 0 + 2b + 4c - 5d = -3 
		 *  
		 * -2x + 5y - z + w + a + 3b + 0 - 2d = 4
		 * 
		 *  x + 0 + z + w + 0 + 2b + c + d = -4
		 */
		
		double [][] system1 = {
				{2,1,-1,-1,1,0,-1,-1},
				{1,0,2,0,-1,-2,2,2},
				{0,-2,5,4,-1,0,3,1},
				{1,1,-7,3,2,1,-1,0},
				{1,1,2,3,-2,2,2,9},
				{0,-3,-2,2,0,2,4,-5},
				{-2,5,-1,1,1,3,0,-2},
				{1,0,1,1,0,2,1,1},
		};
		double[][] system2 = {
				{1},
				{-1},
				{2},
				{-2},
				{3},
				{-3},
				{4},
				{-4},	
		};
		
		Matrix s1 = new Matrix(system1);
		Matrix s2 = new Matrix(system2);
		Matrix temp = s1.gaussJordanElimination(s2);
		return temp;
		
		
	}
	
	public void printNumberNine() {
		double [][] system1 = {
				{2,1,-1,-1,1,0,-1,-1},
				{1,0,2,0,-1,-2,2,2},
				{0,-2,5,4,-1,0,3,1},
				{1,1,-7,3,2,1,-1,0},
				{1,1,2,3,-2,2,2,9},
				{0,-3,-2,2,0,2,4,-5},
				{-2,5,-1,1,1,3,0,-2},
				{1,0,1,1,0,2,1,1},
		};
		
		Matrix A = new Matrix(system1);
		
		double aDeterminant = A.findDeterminant();
		System.out.println("A determinant = " + aDeterminant + "\n");
		
		Matrix aInverse = A.inverse();
		System.out.println("A inverse = \n" + aInverse);
		
		double aInverseDet = aInverse.findDeterminant(); 
		System.out.println("A inverse determinant  " + aInverseDet);
		
		
		double product = aDeterminant * aInverseDet; 
		System.out.println("\nProduct = " + product);
		
		
		
		
		
		
		

		
		
		
		
		

		
		
		
		

		
		
		
		

    } 
	
	public void printConditionNumber() {
		double [][] system1 = {
				{2,1,-1,-1,1,0,-1,-1},
				{1,0,2,0,-1,-2,2,2},
				{0,-2,5,4,-1,0,3,1},
				{1,1,-7,3,2,1,-1,0},
				{1,1,2,3,-2,2,2,9},
				{0,-3,-2,2,0,2,4,-5},
				{-2,5,-1,1,1,3,0,-2},
				{1,0,1,1,0,2,1,1},
		};
		
		Matrix A = new Matrix(system1);
		Matrix AInverse = A.inverse();
		
		double[][] aArray = new double[A.getRows()][A.getColumns()];
		double[][] aInverseArray = new double[A.getRows()][A.getColumns()];
		
		
		////////////////////////////////////////////////////////////// absolute values 
		for(int i = 0; i < A.getRows(); i++) {
			for(int j = 0; j < A.getColumns(); j++) {
				aArray[i][j] = Math.abs(A.getValue(i, j));
			}
		}
		
		for(int i = 0; i < A.getRows(); i++) {
			for(int j = 0; j < A.getColumns(); j++) {
				aInverseArray[i][j] = Math.abs(AInverse.getValue(i, j));
			}
		}
		////////////////////////////////////////////////////////////// absolute values 

		
		////////////////////////////////////////////////////////////// adding all the values in the rows 

		double[] temp = new double[aArray[0].length];
		double[] temp2 = new double[aArray[0].length];

		
		for(int i = 0; i < aArray.length; i++) {
			double rowValue = 0.0;
			for(int j = 0; j < aArray[0].length; j++) {
				rowValue += aArray[i][j];
			}
			temp[i] = rowValue; 
		}
		
		for(int i = 0; i < aInverseArray.length; i++) {
			double rowValue = 0.0;
			for(int j = 0; j < aInverseArray[0].length; j++) {
				rowValue += aInverseArray[i][j];
			}
			temp2[i] = rowValue; 
		}
		
		////////////////////////////////////////////////////////////// adding all the values in the rows 

		
		////////////////////////////////////////////////////////////// finding max 

		
		double condition1 = 0.0; 
		for(int i = 0; i < temp.length; i++) {
			if(condition1 < temp[i]) {
				condition1 = temp[i];
			}
		}
		
		double condition2 = 0.0; 
		for(int i = 0; i < temp2.length; i++) {
			if(condition2 < temp2[i]) {
				condition2 = temp2[i];
			}
		}
		
		System.out.println("Condition 1 = " + condition1);
		System.out.println("Condition 2 = " + condition2);
		System.out.println(condition1 * condition2);
		
	}

}
			

		
		
		
		
