package AlgMatrixProject1;

import java.util.Random;
import java.util.Scanner;

//Taylor Vandenberg
//Design Analysis of Alg
//Proj 1
//Matrix mult 3 ways

public class matrixMult{
    //print a matrix
    static void matrixResult(double A[][], int n){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++)
                System.out.print(A[i][j] + " ");
            System.out.println();
        }
        System.out.println();
    }
    //creates a random array
    static double[][] randomArray(int n){
        double[][] newMatrix = new double [n][n];
        Random r = new Random();
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                Integer randomNum = r.nextInt()% 10;
                newMatrix[i][j] = randomNum;
            }
        }
        return newMatrix;
    }

    //strassen stuff
    //strassen split
    public void matSplit(double[][] matrix1, double[][] matrix2, int A, int B){
        for(int i=0, j=A; i<matrix2.length; i++, j++){
            for(int k=0, l=B; k<matrix2.length; k++, l++){
                matrix2[i][k] = matrix1[j][l];
            }
        }
    }
    //strassen join
    public void matJoin(double[][] matrix2, double[][] matrix1, int A, int B){
        for(int i=0, j=A; i<matrix2.length; i++, j++){
            for(int k=0, l=B; k<matrix2.length; k++, l++){
                matrix1[j][l] = matrix2[i][k];
            }
        }
    }
    //strassen add
    public double[][] matAdd(double[][] matrix1, double[][] matrix2){
        int n = matrix1.length;
        double[][] matrixFinal = new double[n][n];

        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrixFinal[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return matrixFinal;
    }
    //strassen subtract
    public double[][] matSubtract(double[][] matrix1, double[][] matrix2){
        int n=matrix1.length;
        double[][] matrixFinal = new double[n][n];

        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrixFinal[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }
        return matrixFinal;
    }
    //naive div and conquer
    public double[][] naiveConquer(double[][] matrix1, double[][] matrix2){
        int n = matrix1.length;
        double[][] matrixFinal = new double[n][n];
        if(n==1){
            matrixFinal[0][0] = matrix1[0][0] * matrix2[0][0];
        }
        else{
            //matrix devision prep
            double[][] A11 = new double[n/2][n/2];
            double[][] A12 = new double[n/2][n/2];
            double[][] A21 = new double[n/2][n/2];
            double[][] A22 = new double[n/2][n/2];

            double[][] B11 = new double[n/2][n/2];
            double[][] B12 = new double[n/2][n/2];
            double[][] B21 = new double[n/2][n/2];
            double[][] B22 = new double[n/2][n/2];

            //divide matrix1 and 2 into 4 by using the split function
            matSplit(matrix1, A11, 0, 0);
            matSplit(matrix1, A12, 0, n/2);
            matSplit(matrix1, A21, n/2, 0);
            matSplit(matrix1, A22, n/2, n/2);

            matSplit(matrix2, B11, 0, 0);
            matSplit(matrix2, B12, 0, n/2);
            matSplit(matrix2, B21, n/2, 0);
            matSplit(matrix2, B22, n/2, n/2);

            //the maths
            double[][] C11 = matAdd(strassenMult(A11, B11), strassenMult(A12, B21));
            double[][] C12 = matAdd(strassenMult(A11, B12), strassenMult(A12, B22));
            double[][] C21 = matAdd(strassenMult(A21, B11), strassenMult(A22, B21));
            double[][] C22 = matAdd(strassenMult(A21, B12), strassenMult(A22, B22));

            //join back together
            matJoin(C11, matrixFinal, 0, 0);
            matJoin(C12, matrixFinal, 0, n/2);
            matJoin(C21, matrixFinal, n/2, 0);
            matJoin(C22, matrixFinal, n/2, n/2);
        }
        return matrixFinal;
    }
    //strassen main
    public double[][] strassenMult(double[][] matrix1, double[][] matrix2){
        int n = matrix1.length;
        double[][] matrixFinal = new double[n][n];
        //if this part is true no more division needed
        if(n==1){
            matrixFinal[0][0] = matrix1[0][0] * matrix2[0][0];
        }
        //majority of the shenanigans
        else{
            //matrix devision prep
            double[][] A11 = new double[n/2][n/2];
            double[][] A12 = new double[n/2][n/2];
            double[][] A21 = new double[n/2][n/2];
            double[][] A22 = new double[n/2][n/2];

            double[][] B11 = new double[n/2][n/2];
            double[][] B12 = new double[n/2][n/2];
            double[][] B21 = new double[n/2][n/2];
            double[][] B22 = new double[n/2][n/2];

            //divide matrix1 and 2 into 4 by using the split function
            matSplit(matrix1, A11, 0, 0);
            matSplit(matrix1, A12, 0, n/2);
            matSplit(matrix1, A21, n/2, 0);
            matSplit(matrix1, A22, n/2, n/2);

            matSplit(matrix2, B11, 0, 0);
            matSplit(matrix2, B12, 0, n/2);
            matSplit(matrix2, B21, n/2, 0);
            matSplit(matrix2, B22, n/2, n/2);

            //7 mults
            double[][] P1 = strassenMult(A11, (matSubtract(B12, B22)));
            double[][] P2 = strassenMult(matAdd(A11, A12), B22);
            double[][] P3 = strassenMult(matAdd(A21, A22), B11);
            double[][] P4 = strassenMult(A22, matSubtract(B21, B11));
            double[][] P5 = strassenMult(matAdd(A11, A22), matAdd(B11, B22));
            double[][] P6 = strassenMult(matSubtract(A12, A22), matAdd(B21, B22));
            double[][] P7 = strassenMult(matSubtract(A11, A21), matAdd(B11, B12));

            double[][] C11 = matAdd(matAdd(P4, P5), matSubtract(P6, P2));
            double[][] C12 = matAdd(P1, P2);
            double[][] C21 = matAdd(P3, P4);
            double[][] C22 = matAdd(matSubtract(P1, P3), matSubtract(P5, P7));

            //join back together
            matJoin(C11, matrixFinal, 0, 0);
            matJoin(C12, matrixFinal, 0, n/2);
            matJoin(C21, matrixFinal, n/2, 0);
            matJoin(C22, matrixFinal, n/2, n/2);
        }
        return matrixFinal;
    }
    public static void main(String[] args){
        matrixMult staticFix = new matrixMult();
        //n setup
        System.out.println("Input N (2^n): ");
        Scanner kb = new Scanner(System.in);
        int num1 = kb.nextInt();
        System.out.println();
        int n = (int) Math.round(Math.pow(2, num1));

        //random matrix and empty finals created
        double[][] matrix1 = randomArray(n);
        double[][] matrix2 = randomArray(n);
        double[][] bruteMatrix = new double[n][n];
        double[][] naiveMatrixFinal = new double[n][n];
        double[][] strassMatrix = new double[n][n];

        System.out.println("Input Matrices:\n");
        matrixResult(matrix1, n);
        matrixResult(matrix2, n);
        
        //brute force
        long nanoTime1 = System.nanoTime();
        System.out.println("Classical Matrix:\n");
        for(int i=0; i<n; i++){ //loop1
            for(int j=0; j<n; j++){ //loop2
                for (int k=0; k<n; k++){ //loop3
                    bruteMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }//total of n^3
        long nanoTime2 = System.nanoTime();
        long finalNanoTime = nanoTime2 - nanoTime1;
        System.out.println("Time taken: "+ finalNanoTime);
        matrixResult(bruteMatrix, n);

        //naive divide and conquer
        System.out.println("Naive Divide-And-Conquer:\n");
        nanoTime1 = System.nanoTime();
        naiveMatrixFinal = staticFix.naiveConquer(matrix1, matrix2);
        nanoTime2 = System.nanoTime();
        finalNanoTime = nanoTime2 - nanoTime1;
        System.out.println("Time taken: "+ finalNanoTime);
        matrixResult(naiveMatrixFinal, n);

        //Strassen, recursive stuff was done as methods, see above
        System.out.println("Strassen Matrix:\n");
        nanoTime1 = System.nanoTime();
        strassMatrix = staticFix.strassenMult(matrix1, matrix2);
        nanoTime2 = System.nanoTime();
        finalNanoTime = nanoTime2 - nanoTime1;
        System.out.println("Time taken: "+ finalNanoTime);
        matrixResult(strassMatrix, n);


        //sanity check
        double[][] sanityMat1 = {{2,0,-1,6},{3,7,8,0},{-5,1,6,-2},{8,0,2,7}};
        double[][] sanityMat2 = {{0,1,6,3},{-2,8,7,1},{2,0,-1,0},{9,1,6,-2}};
        double[][] sanityMat3 = new double[4][4];
        n = sanityMat1.length;
        System.out.println("\n\n\nSanity Check Matrix\nClassical: \n");
        nanoTime1 = System.nanoTime();
        for(int i=0; i<n; i++){ //loop1
            for(int j=0; j<n; j++){ //loop2
                for (int k=0; k<n; k++){ //loop3
                    sanityMat3[i][j] += sanityMat1[i][k] * sanityMat2[k][j];
                }
            }
        }//total of n^3
        nanoTime2 = System.nanoTime();
        finalNanoTime = nanoTime2 - nanoTime1;
        System.out.println("Time taken: "+ finalNanoTime);
        matrixResult(sanityMat3, n);

        System.out.println("Div & Conquer: \n");
        nanoTime1 = System.nanoTime();
        sanityMat3 = staticFix.naiveConquer(sanityMat1, sanityMat2);
        nanoTime2 = System.nanoTime();
        finalNanoTime = nanoTime2 - nanoTime1;
        System.out.println("Time taken: "+ finalNanoTime);
        matrixResult(sanityMat3, 4);

        System.out.println("Strassen \n");
        nanoTime1 = System.nanoTime();
        sanityMat3 = staticFix.strassenMult(sanityMat1, sanityMat2);
        nanoTime2 = System.nanoTime();
        finalNanoTime = nanoTime2 - nanoTime1;
        System.out.println("Time taken: "+ finalNanoTime);
        matrixResult(sanityMat3, 4);
    }
}