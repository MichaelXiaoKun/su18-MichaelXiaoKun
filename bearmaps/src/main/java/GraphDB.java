import org.xml.sax.SAXException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;


/**
 * Graph for storing all of the intersection (vertex) and road (edge) information.
 * Uses your GraphBuildingHandler to convert the XML files into a graph. Your
 * code must include the vertices, adjacent, distance, closest, lat, and lon
 * methods. You'll also need to include instance variables and methods for
 * modifying the graph (e.g. addNode and addEdge).
 *
 * @author Kevin Lowe, Antares Chen, Kevin Lin
 */
public class GraphDB {
    private Map<Long, Node> vertex = new HashMap<>();
    private Map<Long, Edge> edges = new HashMap<>();
    Map<Long, LocationParams> locations = new HashMap<>();
    private KDTree tree;
    private double sd = Double.MAX_VALUE;
    private long result;

    public class Node {
        long id;
        double longitude;
        double latitude;
        ArrayList<Long> adjacentList = new ArrayList<>();

        Node(long id, double longitude, double latitude) {
            this.id = id;
            this.longitude = longitude;
            this.latitude = latitude;
        }
        public void connectTo(long vertexId) {
            adjacentList.add(vertexId);
        }
    }

    static class Edge {
        private long edgeid;
        private List<Long> vertices;

        Edge(long edgeid, List vertices) {
            this.edgeid = edgeid;
            this.vertices = vertices;
        }
    }


    /**
     * This constructor creates and starts an XML parser, cleans the nodes, and prepares the
     * data structures for processing. Modify this constructor to initialize your data structures.
     * @param dbPath Path to the XML file to be parsed.
     */
    public GraphDB(String dbPath) {
        List<Long> vList = new ArrayList<>();
        File inputFile = new File(dbPath);
        try (FileInputStream inputStream = new FileInputStream(inputFile)) {
            SAXParserFactory factory = SAXParserFactory.newInstance();
            SAXParser saxParser = factory.newSAXParser();
            saxParser.parse(inputStream, new GraphBuildingHandler(this));
        } catch (ParserConfigurationException | SAXException | IOException e) {
            e.printStackTrace();
        }
        clean();
        for (Node t : vertex.values()) {
            vList.add(t.id);
        }
        tree = new KDTree();
        tree.root = tree.kdTreeConstructor(vList, 0);

    }

    /**
     * Helper to process strings into their "cleaned" form, ignoring punctuation and capitalization.
     * @param s Input string.
     * @return Cleaned string.
     */
    private static String cleanString(String s) {
        return s.replaceAll("[^a-zA-Z ]", "").toLowerCase();
    }

    /**
     * Remove nodes with no connections from the graph.
     * While this does not guarantee that any two nodes in the remaining graph are connected,
     * we can reasonably assume this since typically roads are connected.
     */
    private void clean() {
        Set<Long> removedID = new HashSet<>();
        for (Long i : vertex.keySet()) {
            if (vertex.get(i).adjacentList.size() == 0) {
                removedID.add(i);
            }
        }
        for (Long j : removedID) {
            vertex.remove(j);
        }
    }
    public void addEdge(Long id, List<Long> vertices) {
        Edge e = new Edge(id, vertices);
        edges.put(id, e);
        for (int i = 0; i < vertices.size() - 1; i += 1) {
            long firstVertexID = vertices.get(i);
            long secondVertexID = vertices.get(i + 1);
            vertex.get(firstVertexID).connectTo(secondVertexID);
            vertex.get(secondVertexID).connectTo(firstVertexID);
        }
    }

    public void addNode(long id, double longitude, double latitude) {
        vertex.put(id, new Node(id, longitude, latitude));
    }


    /**
     * Returns the longitude of vertex <code>v</code>.
     * @param v The ID of a vertex in the graph.
     * @return The longitude of that vertex, or 0.0 if the vertex is not in the graph.
     */
    double lon(long v) {
        return vertex.get(v).longitude;
    }

    /**
     * Returns the latitude of vertex <code>v</code>.
     * @param v The ID of a vertex in the graph.
     * @return The latitude of that vertex, or 0.0 if the vertex is not in the graph.
     */
    double lat(long v) {
        return vertex.get(v).latitude;
    }

    /**
     * Returns an iterable of all vertex IDs in the graph.
     * @return An iterable of all vertex IDs in the graph.
     */
    Iterable<Long> vertices() {
        return vertex.keySet();
    }

    /**
     * Returns an iterable over the IDs of all vertices adjacent to <code>v</code>.
     * @param v The ID for any vertex in the graph.
     * @return An iterable over the IDs of all vertices adjacent to <code>v</code>, or an empty
     * iterable if the vertex is not in the graph.
     */
    Iterable<Long> adjacent(long v) {
        return vertex.get(v).adjacentList;
    }

    /**
     * Returns the great-circle distance between two vertices, v and w, in miles.
     * Assumes the lon/lat methods are implemented properly.
     * @param v The ID for the first vertex.
     * @param w The ID for the second vertex.
     * @return The great-circle distance between vertices and w.
     * @source https://www.movable-type.co.uk/scripts/latlong.html
     */
    public double distance(long v, long w) {
        double phi1 = Math.toRadians(lat(v));
        double phi2 = Math.toRadians(lat(w));
        double dphi = Math.toRadians(lat(w) - lat(v));
        double dlambda = Math.toRadians(lon(w) - lon(v));

        double a = Math.sin(dphi / 2.0) * Math.sin(dphi / 2.0);
        a += Math.cos(phi1) * Math.cos(phi2) * Math.sin(dlambda / 2.0) * Math.sin(dlambda / 2.0);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return R * c;
    }

    /**
     * Returns the ID of the vertex closest to the given longitude and latitude.
     * @param lon The given longitude.
     * @param lat The given latitude.
     * @return The ID for the vertex closest to the <code>lon</code> and <code>lat</code>.
     */
    private static double euclidean(double x1, double x2, double y1, double y2) {
        return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
    }


    public class KDTree {
        private TreeNode root;
        private class TreeNode {
            private long id;
            private TreeNode left;
            private TreeNode right;
            private int depth;

            TreeNode(long id, int depth) {
                this.id = id;
                this.left = null;
                this.right = null;
                this.depth = depth;
            }

            public void helper(double lon, double lat) {
                KDTree.TreeNode tmp = this;
                double td = euclidean(projectToX(lon(tmp.id), lat(tmp.id)), projectToX(lon, lat),
                                      projectToY(lon(tmp.id), lat(tmp.id)), projectToY(lon, lat));
                if (td < sd) {
                    sd = td;
                    result = tmp.id;
                }
                if (depth % 2 == 0) {
                    if (projectToX(lon, lat) < projectToX(lon(tmp.id), lat(tmp.id))) {
                        if (tmp.left != null) {
                            left.helper(lon, lat);
                        }
                        if (sd > Math.abs(projectToX(lon, lat) - projectToX(lon(id), lat(id)))) {
                            if (tmp.right != null) {
                                right.helper(lon, lat);
                            }
                        }
                    } else {
                        if (tmp.right != null) {
                            right.helper(lon, lat);
                        }
                        if (sd > Math.abs(projectToX(lon, lat) - projectToX(lon(id), lat(id)))) {
                            if (tmp.left != null) {
                                left.helper(lon, lat);
                            }
                        }
                    }
                } else {
                    if (projectToY(lon, lat) < projectToY(lon(tmp.id), lat(tmp.id))) {
                        if (tmp.left != null) {
                            left.helper(lon, lat);
                        }
                        if (sd > Math.abs(projectToY(lon, lat) - projectToY(lon(id), lat(id)))) {
                            if (tmp.right != null) {
                                right.helper(lon, lat);
                            }
                        }
                    } else {
                        if (tmp.right != null) {
                            right.helper(lon, lat);
                        }
                        if (sd > Math.abs(projectToY(lon, lat) - projectToY(lon(id), lat(id)))) {
                            if (tmp.left != null) {
                                left.helper(lon, lat);
                            }
                        }
                    }
                }
            }
        }

        public TreeNode kdTreeConstructor(List<Long> listsOfPoints, int depth) {
            if (listsOfPoints.isEmpty()) {
                return null;
            }

            int axis = depth % 2;
            if (axis == 0) {
                Collections.sort(listsOfPoints, (u, v) -> Double.compare(projectToX(lon(u), lat(u)),
                        projectToX(lon(v), lat(v))));
            }
            if (axis == 1) {
                Collections.sort(listsOfPoints, (u, v) -> Double.compare(projectToY(lon(u), lat(u)),
                        projectToY(lon(v), lat(v))));
            }

            TreeNode trees = new TreeNode(listsOfPoints.get(listsOfPoints.size() / 2), depth);
            trees.left = kdTreeConstructor(listsOfPoints.subList(0, listsOfPoints.size() / 2),
                    depth + 1);
            trees.right = kdTreeConstructor(listsOfPoints.subList(listsOfPoints.size() / 2 + 1,
                    listsOfPoints.size()), depth + 1);
            return trees;
        }
    }

    public long closest(double lon, double lat) {
        sd = euclidean(projectToX(lon(tree.root.id), lat(tree.root.id)),
                projectToX(lon, lat), projectToY(lon(tree.root.id),
                        lat(tree.root.id)), projectToY(lon, lat));
        result = tree.root.id;
        sd = euclidean(projectToX(lon(tree.root.id), lat(tree.root.id)),
                projectToX(lon, lat), projectToY(lon(tree.root.id),
                        lat(tree.root.id)), projectToY(lon, lat));
        tree.root.helper(lon, lat);
        return result;
    }


    /**
     * Return the Euclidean x-value for some point, p, in Berkeley. Found by computing the
     * Transverse Mercator projection centered at Berkeley.
     * @param lon The longitude for p.
     * @param lat The latitude for p.
     * @return The flattened, Euclidean x-value for p.
     * @source https://en.wikipedia.org/wiki/Transverse_Mercator_projection
     */

    static double projectToX(double lon, double lat) {
        double dlon = Math.toRadians(lon - ROOT_LON);
        double phi = Math.toRadians(lat);
        double b = Math.sin(dlon) * Math.cos(phi);
        return (K0 / 2) * Math.log((1 + b) / (1 - b));
    }

    /**
     * Return the Euclidean y-value for some point, p, in Berkeley. Found by computing the
     * Transverse Mercator projection centered at Berkeley.
     * @param lon The longitude for p.
     * @param lat The latitude for p.
     * @return The flattened, Euclidean y-value for p.
     * @source https://en.wikipedia.org/wiki/Transverse_Mercator_projection
     */
    static double projectToY(double lon, double lat) {
        double dlon = Math.toRadians(lon - ROOT_LON);
        double phi = Math.toRadians(lat);
        double con = Math.atan(Math.tan(phi) / Math.cos(dlon));
        return K0 * (con - Math.toRadians(ROOT_LAT));
    }

    /**
     * In linear time, collect all the names of OSM locations that prefix-match the query string.
     * @param prefix Prefix string to be searched for. Could be any case, with our without
     *               punctuation.
     * @return A <code>List</code> of the full names of locations whose cleaned name matches the
     * cleaned <code>prefix</code>.
     */
    public List<String> getLocationsByPrefix(String prefix) {
        return Collections.emptyList();
    }

    /**
     * Collect all locations that match a cleaned <code>locationName</code>, and return
     * information about each node that matches.
     * @param locationName A full name of a location searched for.
     * @return A <code>List</code> of <code>LocationParams</code> whose cleaned name matches the
     * cleaned <code>locationName</code>
     */
    public List<LocationParams> getLocations(String locationName) {
        return Collections.emptyList();
    }

    /**
     * Returns the initial bearing between vertices <code>v</code> and <code>w</code> in degrees.
     * The initial bearing is the angle that, if followed in a straight line along a great-circle
     * arc from the starting point, would take you to the end point.
     * Assumes the lon/lat methods are implemented properly.
     * @param v The ID for the first vertex.
     * @param w The ID for the second vertex.
     * @return The bearing between <code>v</code> and <code>w</code> in degrees.
     * @source https://www.movable-type.co.uk/scripts/latlong.html
     */
    double bearing(long v, long w) {
        double phi1 = Math.toRadians(lat(v));
        double phi2 = Math.toRadians(lat(w));
        double lambda1 = Math.toRadians(lon(v));
        double lambda2 = Math.toRadians(lon(w));

        double y = Math.sin(lambda2 - lambda1) * Math.cos(phi2);
        double x = Math.cos(phi1) * Math.sin(phi2);
        x -= Math.sin(phi1) * Math.cos(phi2) * Math.cos(lambda2 - lambda1);
        return Math.toDegrees(Math.atan2(y, x));
    }


    /** Radius of the Earth in miles. */
    private static final int R = 3963;
    /** Latitude centered on Berkeley. */
    private static final double ROOT_LAT = (MapServer.ROOT_ULLAT + MapServer.ROOT_LRLAT) / 2;
    /** Longitude centered on Berkeley. */
    private static final double ROOT_LON = (MapServer.ROOT_ULLON + MapServer.ROOT_LRLON) / 2;
    /**
     * Scale factor at the natural origin, Berkeley. Prefer to use 1 instead of 0.9996 as in UTM.
     * @source https://gis.stackexchange.com/a/7298
     */
    private static final double K0 = 1.0;

//    public GraphDB() {
//        List<Long> vlist = new ArrayList<>(vertex.keySet());
//        tree = new KDTree(vlist);
//    }

//    public static void main(String[] args) {
//        GraphDB g = new GraphDB();
//        g.addNode(1L, 27.4848, 5.2345);
//        g.addNode(23, 27.1, 5.1);
//        g.addNode(2L, 30.3444444, 4.93456);
//        g.addNode(3L, 23.456543, 3.394543);
//        System.out.println(g.closest(27, 5));
//    }
}
