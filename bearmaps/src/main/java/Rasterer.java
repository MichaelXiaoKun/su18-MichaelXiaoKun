/**
 * This class provides all code necessary to take a query box and produce
 * a query result. The getMapRaster method must return a Map containing all
 * seven of the required fields, otherwise the front end code will probably
 * not draw the output correctly.
 */
public class Rasterer {
    /** The max image depth level. */
    public static final int MAX_DEPTH = 7;
    /**
     * Takes a user query and finds the grid of images that best matches the query. These images
     * will be combined into one big image (rastered) by the front end. The grid of images must obey
     * the following properties, where image in the grid is referred to as a "tile".
     * <ul>
     *     <li>The tiles collected must cover the most longitudinal distance per pixel (LonDPP)
     *     possible, while still covering less than or equal to the amount of longitudinal distance
     *     per pixel in the query box for the user viewport size.</li>
     *     <li>Contains all tiles that intersect the query bounding box that fulfill the above
     *     condition.</li>
     *     <li>The tiles must be arranged in-order to reconstruct the full image.</li>
     * </ul>
     * @param params The RasterRequestParams containing coordinates of the query box and the browser
     *               viewport width and height.
     * @return A valid RasterResultParams containing the computed results.
     */
    public RasterResultParams getMapRaster(RasterRequestParams params) {
        double lonDPP = lonDPP(params.lrlon, params.ullon, params.w);
        if (params.lrlon < MapServer.ROOT_ULLON || params.lrlat > MapServer.ROOT_ULLAT) {
            System.out.println("The latitude and longitude are wrong!");
            return RasterResultParams.queryFailed();
        }
        int depth = 0;
        for (int i = 7; i >= 0; i--) {
            if (lonDPP < MapServer.ROOT_LONDPP / Math.pow(2, i - 1)) {
                depth = i;
                break;
            }
        }
        int markULLongitude = 0, markULLatitude = 0;
        int markLRLongitude = (int) Math.pow(2, depth) - 1;
        int markLRLatitude = (int) Math.pow(2, depth) - 1;
        int gridSize1f = 0, gridSize2f = 0;
        for (int i = 0; i < Math.pow(2, depth); i++) {
            double ullongitude = MapServer.ROOT_ULLON
                    + (i * MapServer.ROOT_LON_DELTA) / Math.pow(2, depth);
            if (ullongitude > params.ullon) {
                markULLongitude = Math.max(markULLongitude, i - 1);
                break;
            }
        }
        for (int j = 0; j < Math.pow(2, depth); j++) {
            double ullatitude = MapServer.ROOT_ULLAT
                    - (j * MapServer.ROOT_LAT_DELTA) / Math.pow(2, depth);
            if (ullatitude < params.ullat) {
                markULLatitude = Math.max(markULLatitude, j - 1);
                break;
            }
        }
        for (int j = 0; j < Math.pow(2, depth); j++) {
            double lrlatitude = MapServer.ROOT_LRLAT
                    + (j * MapServer.ROOT_LAT_DELTA) / Math.pow(2, depth);
            if (lrlatitude > params.lrlat) {
                gridSize1f = Math.min(markLRLatitude, (int) Math.pow(2, depth) - j);
                markLRLatitude = Math.max(0, j - 1);
                break;
            }
        }
        for (int i = 0; i < Math.pow(2, depth); i++) {
            double lrlongitude = MapServer.ROOT_LRLON
                    - (i * MapServer.ROOT_LON_DELTA) / Math.pow(2, depth);
            if (lrlongitude < params.lrlon) {
                gridSize2f = Math.min(markLRLongitude, (int) Math.pow(2, depth) - i);
                markLRLongitude = Math.max(0, i - 1);
                break;
            }
        }
        int a = Math.abs(gridSize1f - markULLatitude + 1);
        int b = Math.abs(gridSize2f - markULLongitude + 1);
        String[][] grid = new String[a][b];
        double ullongitude = MapServer.ROOT_ULLON
                + (markULLongitude * (MapServer.ROOT_LON_DELTA)) / Math.pow(2, depth);
        double ullatitude = MapServer.ROOT_ULLAT
                - (markULLatitude * (MapServer.ROOT_LAT_DELTA)) / Math.pow(2, depth);
        double lrlongitude = MapServer.ROOT_LRLON
                - ((markLRLongitude * (MapServer.ROOT_LON_DELTA))) / Math.pow(2, depth);
        double lrlatitude = MapServer.ROOT_LRLAT
                + (markLRLatitude * (MapServer.ROOT_LAT_DELTA)) / Math.pow(2, depth);
        int countX = markULLongitude;
        int countY = markULLatitude;
        for (int i = 0; i < grid.length; i++) {
            for (int j = 0; j < grid[0].length; j++) {
                grid[i][j] = "d" + depth + "_x" + (j + countX) + "_y" + (countY + i) + ".png";
            }
        }
        RasterResultParams.Builder result = new RasterResultParams.Builder();
        result.setRasterLrLat(lrlatitude);
        result.setRasterLrLon(lrlongitude);
        result.setRasterUlLat(ullatitude);
        result.setRasterUlLon(ullongitude);
        result.setDepth(depth);
        result.setRenderGrid(grid);
        result.setQuerySuccess(true);
        return result.create();
    }

    /**
     * Calculates the lonDPP of an image or query box
     * @param lrlon Lower right longitudinal value of the image or query box
     * @param ullon Upper left longitudinal value of the image or query box
     * @param width Width of the query box or image
     * @return lonDPP
     */
    private double lonDPP(double lrlon, double ullon, double width) {
        return (lrlon - ullon) / width;
    }

    public static void main(String[] args) {
        RasterRequestParams.Builder a = new RasterRequestParams.Builder();
        Rasterer b;
        b = new Rasterer();
        a.setLrlat(37.87548268822065);
        a.setLrlon(-122.24053369025242);
        a.setUllat(37.87655856892288);
        a.setUllon(-122.24163047377972);
        a.setH(875.0);
        a.setW(892.0);
        RasterResultParams c = b.getMapRaster(a.create());
        System.out.println(c.rasterLrLat);
        System.out.println(c.rasterLrLon);
        System.out.println(c.rasterUlLat);
        System.out.println(c.rasterUlLon);


        for (int i = 0; i < c.renderGrid.length; i++) {
            for (int j = 0; j < c.renderGrid[0].length; j++) {
                System.out.println(c.renderGrid[i][j]);
            }
        }


        System.out.println(c.renderGrid.length);
        System.out.println(c.renderGrid[0].length);


    }
}
