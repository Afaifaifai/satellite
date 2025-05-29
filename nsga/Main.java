public class Main {
    public static void main(String[] a) {
        Settings s = new Settings();
        int sequence_num = 1;
        String[][] s_access_states = Csv2Array.readCsvToArray(String.format("%s%s_%03d.%s", s.dir, s.iridium_access_states, sequence_num, s.file_suffix));
        String[][] s_coords = Csv2Array.readCsvToArray(String.format("%s%s_%03d.%s", s.dir, s.iridium_coords, sequence_num, s.file_suffix));
        String[][] s_distances = Csv2Array.readCsvToArray(String.format("%s%s_%03d.%s", s.dir, s.iridium_distances, sequence_num, s.file_suffix));

        double[][] coords = Csv2Array.to_double_array(s_coords), distances = Csv2Array.to_double_array(s_distances);
        int[][] access_states = Csv2Array.to_int_array(s_access_states);
        // show(access_states);

        

    }

    public static void show(double[][] ary) {
        for (double[] temp : ary) {
            for (double tt : temp)
                System.out.printf("%.2f ", tt);
            System.out.println();
        }
        System.out.println();
    }

    static void show(int[][] ary) {
        for (int[] temp : ary) {
            for (int tt : temp)
                System.out.printf("%d ", tt);
            System.out.println();
        }
        System.out.println();
    }
}
