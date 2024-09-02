import java.util.Arrays;
import java.util.Scanner;

public class Matrix {

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        System.out.print("Enter no of rows or columns of a square matrix: ");
        int n = sc.nextInt();
        int[][] A = new int[n][n];
        int[][] B = new int[n][n];
        System.out.println("Enter all "+(n*n)+" elements in A matrix:");
        inputMatrix(A);
        System.out.println("Enter all "+(n*n)+" elements in B matrix:");
        inputMatrix(B);

        // int[][] A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        // int[][] B = {{3, 5, 6}, {2, 5, 9}, {7, 8, 1}};

        System.out.println("\nMatrix A:");
        toPrint(A);
        System.out.println("\nMatrix B:");
        toPrint(B);
        System.out.println("\nMultiplication of A and B is:");
        toPrint(toMultiply(A, B));
        System.out.println("\nDifference of A and B is:");
        toPrint(toSubtract(A, B));
    }

    static int[][] inputMatrix(int[][] arr){
        Scanner sc = new Scanner(System.in);
        for (int[] arr1 : arr) {
            for (int j = 0; j<arr.length; j++) {
                arr1[j] = sc.nextInt();
            }
        }
        return arr;
    }

    static void toPrint(int[][] arr) {
        for (int[] arr1 : arr) {
            System.out.print(Arrays.toString(arr1));
            System.out.println("");
        }
    }

    static int[][] toMultiply(int[][] arr1,int[][] arr2){
        int[][] result = new int[arr1.length][arr1.length];
        for(int i=0;i<arr1.length;i++){
            for(int j=0;j<arr1.length;j++){
                result[i][j] = 0;
                for(int k=0;k<arr1.length;k++){
                    result[i][j] = result[i][j] + (arr1[i][k]*arr2[k][j]);
                }
            }
        }
        return result;
    }

    static int[][] toSubtract(int[][] arr1,int[][] arr2){
        int[][] result = new int[arr1.length][arr1.length];
        for(int i=0;i<arr1.length;i++){
            for(int j=0;j<arr1.length;j++){
                result[i][j] = arr1[i][j] - arr2[i][j];
            }
        }
        return result;
    }
}
