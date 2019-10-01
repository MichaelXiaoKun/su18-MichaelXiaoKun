package gitlet;

import java.io.Serializable;
import java.io.File;

public class Blob implements Serializable {
    private String filename;
    private byte[] content;

    private Blob() {

    }

    public Blob(String fileName, byte[] content) {
        this.filename = fileName;
        this.content = content;
    }

    public static Blob createBlobFromFile(String filename) {
        File file = new File(filename);
        if (!file.exists()) {
            System.out.println("File does not exist.");
            return null;
        }
        Blob blob = new Blob();
        blob.filename = filename;
        blob.content = Utils.readContents(file);
        return blob;
    }

    public String getHashCode() {
        return Utils.sha1(content);
    }

    public String getBlobFilename() {
        return getHashCode() + "." + getFilename();
    }

    public byte[] getContent() {
        return content;
    }

    public String getFilename() {
        return filename;
    }

    public void checkout() {
        Utils.writeContents(new File(filename), content);
    }
}
