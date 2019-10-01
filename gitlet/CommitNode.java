package gitlet;

import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.Date;
//import java.util.*;

public class CommitNode implements Serializable {
    private String message;
    private HashMap<String, String> blobs;
    private long dateTime;
    private String prev;

    public CommitNode(String msg, HashMap<String, String> blobs, CommitNode prev) {
        message = msg;
        this.blobs = blobs;
        this.prev = prev == null ? null : prev.getHashcode();
        dateTime = System.currentTimeMillis();
    }

    public String getHashcode() {
        return Utils.sha1(Long.toString(dateTime), message);
    }

    public Date getDateTime() {
        return new Date(dateTime);
    }

    public String convertDateTime() {
        return new SimpleDateFormat("uuuu-MM-dd HH:mm:ss").format(getDateTime());
    }

    public String getMessage() {
        return message;
    }

    public String getBlobFilename(String filename) {
        return blobs.get(filename);
    }

    public String getPreviousCommitNodeFilename() {
        return prev;
    }

    public boolean containsFile(String filename) {
        return blobs.containsKey(filename);
    }

    public HashMap<String, String> getBlobs() {
        return blobs;
    }

    public String getLog() {
        String id = this.getHashcode();
        String time = this.convertDateTime();
        String message1 = this.getMessage();
        return "===\nCommit " + id + "\n" + time + "\n" + message1 + "\n\n";
    }

}
