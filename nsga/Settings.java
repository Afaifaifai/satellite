import java.io.BufferedReader;   // BufferedReader（緩衝讀取器）
import java.io.FileReader;       // FileReader（檔案讀取器）
import java.io.IOException;
import java.util.ArrayList;      // ArrayList（動態陣列）
import java.util.List;


public class Settings {
    String dir = "../data/";
    String iridium_access_states = "iridium_access_states";
    String iridium_coords = "iridium_coords";
    String iridium_distances = "iridium_distances";
    String file_suffix = "csv";

    // GA
    int num_population = 100;
}


class Csv2Array {
    public static String[][] readCsvToArray(String filePath) {
        List<String[]> rows = new ArrayList<>();  // 用來暫存每一行資料

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                // 假設 CSV 欄位以逗號分隔，並沒有考慮到引號內含逗號的複雜情況
                String[] values = line.split(",");
                rows.add(values);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // 將 List<String[]> 轉成 String[][]
        String[][] array = new String[rows.size()][];
        for (int i = 0; i < rows.size(); i++) {
            array[i] = rows.get(i);
        }
        return array;
    }

    public static int[][] to_int_array(String[][] ary) {
        int len = ary.length;
        int[][] result = new int[len][];
        for (int i = 0; i < len ; i++) {
            int[] temp = new int[ary[i].length];
            for (int j = 0; j < temp.length; j++)
                temp[j] = Integer.parseInt(ary[i][j]);
            result[i] = temp;
        }

        return result;
    }
    public static double[][] to_double_array(String[][] ary) {
        int len = ary.length;
        double[][] result = new double[len][];
        for (int i = 0; i < len ; i++) {
            double[] temp = new double[ary[i].length];
            for (int j = 0; j < temp.length; j++)
                temp[j] = Double.parseDouble(ary[i][j]);
            result[i] = temp;
        }

        return result;
    }
}