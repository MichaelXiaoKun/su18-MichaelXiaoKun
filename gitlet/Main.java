package gitlet;

import java.io.Serializable;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.List;
import java.util.TreeSet;


/* Driver class for Gitlet, the tiny stupid version-control system.
   @author
import java.io.Serializable;
import java.io.File;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.List;
import java.util.TreeSet;
*/
public class Main implements Serializable {
    private CurrentState currentState;
    /* Usage: java gitlet.Main ARGS, where ARGS contains
       <COMMAND> <OPERAND> .... */
    private TreeSet<String> stagedset;
    private TreeSet<String> removedset;

    public CommitNode getCommitNode(String hashCode) {
        if (hashCode == null) {
            return null;
        }
        return (CommitNode) deserialize(hashCode + ".commit");
    }

    public void saveCommitNode(CommitNode node) {
        serialize(node, node.getHashcode() + ".commit");
    }

    public ArrayList<CommitNode> getAllCommitNode() {
        ArrayList<CommitNode> result = new ArrayList<>();
        for (String string : currentState.getCommitList()) {
            result.add(getCommitNode(string));
        }
        return result;
    }

    public void globalLog() {
        ArrayList<CommitNode> commitNodeList = getAllCommitNode();
        for (CommitNode node : commitNodeList) {
            System.out.print(node.getLog());
        }
    }

    public Blob getBlob(String filename) {
        if (filename == null) {
            return null;
        }
        return (Blob) deserialize(filename);
    }

    public void saveBlob(Blob blob) {
        serialize(blob, blob.getBlobFilename());
    }

    public void loadCurrentState() {
        currentState = (CurrentState) deserialize("currentState");
    }

    public void saveCurrentState() {
        serialize(currentState, "currentState");
    }


    public CommitNode getCurrentCommit() {
        if (currentState == null) {
            return null;
        }

        return getCommitNode(currentState.getCurrentCommit());
    }

    public void setCurrentCommit(CommitNode node) {
        if (currentState == null) {
            return;
        }

        currentState.setCurrentCommit(node.getHashcode());
    }

    public static void serialize(Object obj, String path) {
        File file = new File(".gitlet");
        File outFile = new File(file, path + ".ser");
        try {
            ObjectOutputStream out = new ObjectOutputStream(
                    new FileOutputStream(outFile));
            out.writeObject(obj);
            out.close();
        } catch (IOException e) {
            return;
        }
    }

    public static Object deserialize(String path) {
        Object obj = null;
        File file = new File(".gitlet");
        File inFile = new File(file, path + ".ser");
        try {
            ObjectInputStream inp =
                    new ObjectInputStream(new FileInputStream(inFile));
            obj = inp.readObject();
            inp.close();
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
            obj = null;
        }
        return obj;
    }


    public void init() {
        File file = new File(".gitlet");
        if (file.exists()) {
            System.out.println("A gitlet version-control system "
                    + "already exists in the current directory.");
            return;
        }
        file.mkdir();
        currentState = new CurrentState();
        CommitNode initialCommit = new CommitNode("initial commit", new HashMap<>(), null);
        saveCommitNode(initialCommit);
        currentState.setCurrentCommit(initialCommit.getHashcode());
        currentState.putBranch("master", initialCommit.getHashcode());
        currentState.setCurrentBranch("master");
        currentState.putCommitList(initialCommit
                .getHashcode());
        saveCurrentState();
    }

    public void add(String filename) {
        File file = new File(filename);
        if (!file.exists()) {
            System.out.println("File does not exist.");
            return;
        }
        CommitNode currentCommit = getCurrentCommit();
        Blob currentVersionBlob = Blob.createBlobFromFile(filename);
        if (currentCommit.containsFile(filename)) {
            if (currentVersionBlob.getBlobFilename().
                    equals(currentState.getBlobFilename(filename))) {
                // Same Content
                return;
            }
        }
        saveBlob(currentVersionBlob);
        currentState.addBlob(currentVersionBlob.getFilename(),
                currentVersionBlob.getBlobFilename());
    }

    public void commit(String msg) {
        if (msg == null || msg.isEmpty()) {
            System.out.println("Please enter a commit message.");
            return;
        }
        File exist = new File(".gitlet/StagingArea");
        CommitNode currentNode = getCurrentCommit();
        if (currentState.getBlobs().equals(currentNode.getBlobs())) {
            System.out.println("No changes added to the commit.");
            return;
        }
        CommitNode newNode = new CommitNode(msg, currentState.getBlobs(), getCurrentCommit());
        saveCommitNode(newNode);
        currentState.putBranch(currentState.getCurrentBranchTitle(), newNode.getHashcode());
        currentState.setCurrentCommit(newNode.getHashcode());
        currentState.putCommitList(newNode.getHashcode());
    }

    public void remove(String fileName) {
        CommitNode currentCommit = getCurrentCommit();
        if (currentCommit.getBlobFilename(fileName) == null) {
            if (currentState.getBlobFilename(fileName) == null) {
                System.out.println("No reason to remove the file.");
                return;
            }
        } else {
            File f = new File(fileName);
            f.delete();
        }
        currentState.removeBlob(fileName);

    }

    public void log() {
        CommitNode commit = getCurrentCommit();
        while (commit != null) {
            System.out.print(commit.getLog());
            commit = getCommitNode(commit.getPreviousCommitNodeFilename());
        }
    }

    public void reset(String commitID) {
        if (!currentState.getCommitList().contains(commitID)) {
            System.out.println("No commit with that id exists.");
            return;
        }
        CommitNode node = getCommitNode(commitID);
        for (String filename : node.getBlobs().keySet()) {
            if ((!currentState.getBlobs().containsKey(filename))
                    && new File(filename).exists()) {
                System.out.println("There is an untracked file in the way; "
                        + "delete it or add it first.");
                return;
            }
        }
        currentState.getBlobs().clear();
        HashMap<String, String> map = node.getBlobs();
        for (String filename : getCurrentCommit().getBlobs().keySet()) {
            if (!map.containsKey(filename)) {
                File file = new File(filename);
                file.delete();
            }
        }

        for (String filename : map.keySet()) {
            getBlob(map.get(filename)).checkout();
            currentState.getBlobs().put(filename, map.get(filename));
        }
        currentState.setCurrentCommit(node.getHashcode());
        currentState.putBranch(currentState.getCurrentBranchTitle(), commitID);
    }


    public void branch(String branchName) {
        if (currentState.containsBranch(branchName)) {
            System.out.println("A branch with that name already exists.");
        }
        currentState.putBranch(branchName, currentState.getCurrentCommit());
    }

    public void rmBranch(String branchName) {
        if (!currentState.containsBranch(branchName)) {
            System.out.println("A branch with that name does not exist.");
        } else if (currentState.getCurrentBranchTitle().equals(branchName)) {
            System.out.println("Cannot remove the current branch.");
        }
        currentState.deleteBranch(branchName);
    }

    public String getSplitPoint(String branchName1, String branchName2) {
        //get all node on one branch
        String commitNode1 = currentState.getBranchHeadCommitNode(branchName1);
        String commitNode2 = currentState.getBranchHeadCommitNode(branchName2);
        Set<String> commitSet1 = new HashSet<>();
        while (commitNode1 != null) {
            commitSet1.add(commitNode1);
            commitNode1 = getCommitNode(commitNode1).getPreviousCommitNodeFilename();
        }
        while (commitNode2 != null) {
            if (commitSet1.contains(commitNode2)) {
                return commitNode2;
            }
            commitNode2 = getCommitNode(commitNode2).getPreviousCommitNodeFilename();
        }
        return null;
    }

    public void merge(String targetBranchName) {

        if (!currentState.getAllBranchNames().contains(targetBranchName)) {
            System.out.println("A branch with that name does not exist.");
            return;
        }

        if (currentState.getCurrentBranchTitle().equals(targetBranchName)) {
            System.out.println("Cannot merge a branch with itself.");
            return;
        }

        CommitNode splitPoint = getCommitNode(
                getSplitPoint(currentState.getCurrentBranchTitle(),
                        targetBranchName));
        CommitNode targetPoint = getCommitNode(currentState.
                getBranchHeadCommitNode(targetBranchName));
        CommitNode currPoint = getCurrentCommit();
        if (splitPoint.getHashcode().equals(targetPoint.getHashcode())) {
            System.out.println("Given branch is an ancestor of the current branch.");
            return;
        }
        if (splitPoint.getHashcode().equals(currPoint.getHashcode())) {
            currentState.setCurrentCommit(
                    currentState.getBranchHeadCommitNode(targetBranchName));
            currentState.setCurrentBranch(targetBranchName);
            System.out.println("Current branch fast-forwarded.");
            return;
        }
        if (!currPoint.getBlobs().equals(currentState.getBlobs())) {
            System.out.println("You have uncommitted changes.");
            return;
        }

        HashMap<String, String> currentCommitBlobs = currPoint.getBlobs();
        HashMap<String, String> targetBranchBlobs = targetPoint.getBlobs();
        HashMap<String, String> splitPointBlobs = splitPoint.getBlobs();


        HashSet<String> allBlobs = new HashSet<>();
        allBlobs.addAll(currentCommitBlobs.keySet());
        allBlobs.addAll(targetBranchBlobs.keySet());
        allBlobs.addAll(splitPointBlobs.keySet());

        HashMap<String, String> newBlobs = new HashMap<>();
        boolean isConflictHappened = false;

        HashMap<String, String> filesToUpdate = new HashMap<>();
        for (String file : allBlobs) {
            String currentCommitHash = currentCommitBlobs.get(file);
            String targetCommitHash = targetBranchBlobs.get(file);
            String splitPointHash = splitPointBlobs.get(file);
            String resultHash = null;


            boolean currentCommitChanged = isChanged(currentCommitHash, splitPointHash);
            boolean targetCommitChanged = isChanged(targetCommitHash, splitPointHash);
            boolean toStage = true;

            if (currentCommitChanged && targetCommitChanged
                    && isChanged(currentCommitHash, targetCommitHash)) {
                resultHash = mergeConflictCase(currentCommitHash, targetCommitHash);
                isConflictHappened = true;
                toStage = false;
            } else if (currentCommitChanged) {
                resultHash = currentCommitHash;
            } else {
                resultHash = targetCommitHash;
            }

            if (resultHash == null) {
                if (currentCommitHash != null) {
                    filesToUpdate.put(file, null);
                    if (currentCommitHash == null && new File(file).exists()) {
                        System.out.println("There is an untracked file in the way;"
                                + " delete it or add it first.");
                        return;
                    }
                }
            } else {
                if (toStage) {
                    newBlobs.put(file, resultHash);
                } else if (currentCommitBlobs.containsKey(file)) {
                    newBlobs.put(file, currentCommitBlobs.get(file));
                }
                if (!resultHash.equals(currentCommitHash)) {
                    filesToUpdate.put(file, resultHash);
                    if (currentCommitHash == null && new File(file).exists()) {
                        System.out.println("There is an untracked file in the way; "
                                + "delete it or add it first.");
                        return;
                    }
                }
            }
        }

        for (String file : filesToUpdate.keySet()) {
            if (filesToUpdate.get(file) == null) {
                File resultFile = new File(file);
                resultFile.delete();
            } else {
                getBlob(filesToUpdate.get(file)).checkout();
            }
        }
        currentState.setBlobs(newBlobs);
        if (!isConflictHappened) {
            commit("Merged " + currentState.getCurrentBranchTitle()
                    + " with " + targetBranchName + ".");
        } else {
            System.out.println("Encountered a merge conflict.");
        }
    }

    public static boolean isChanged(String current, String split) {
        if (split == null) {
            //System.out.println("SPLIT IS NULL");
            return current != null;
        } else {
            if (current == null) {
                return true;
            }
        }
        return !current.equals(split);
    }

    public static void appendArray(byte[] src, byte[] dest, int destPosition) {
        System.arraycopy(src, 0, dest, destPosition, src.length);
    }

    //when modified in different way

    public String mergeConflictCase(String first, String second) {
        Blob firstBlob = getBlob(first);
        Blob secondBlob = getBlob(second);
        byte[][] contents = {
                "<<<<<<< HEAD\n".getBytes(), firstBlob != null ? firstBlob.getContent()
                : new byte[0],
                "=======\n".getBytes(), secondBlob != null ? secondBlob.getContent() : new byte[0],
                ">>>>>>>\n".getBytes()
        };
        int totalLength = 0;
        for (byte[] content : contents) {
            totalLength += content.length;
        }
        byte[] resultContent = new byte[totalLength];
        int destPosition = 0;
        for (byte[] content : contents) {
            appendArray(content, resultContent, destPosition);
            destPosition += content.length;
        }
        String filename = firstBlob != null
                ? firstBlob.getFilename() : secondBlob.getFilename();
        //System.out.println("Filename: " + filename);
        //System.out.println("File1:\n" + new String(contents[1]));
        //System.out.println("File2:\n" + new String(contents[3]));
        //System.out.println("Result:\n" + new String(resultContent));
        Blob resultBlob = new Blob(filename, resultContent);
        saveBlob(resultBlob);
        return resultBlob.getBlobFilename();
    }

    /////
    public void checkout(String mark, String filename) {
        checkout(currentState.getCurrentCommit(), mark, filename);
    }

    public void checkout(String commitID, String mark, String filename) {
        if (!mark.equals("--")) {
            System.out.println("Incorrect operands.");
            return;
        }
        String completeCommitId = null;

        for (String completeCommitIdCandidate : currentState.getCommitList()) {
            if (completeCommitIdCandidate.startsWith(commitID)) {
                if (completeCommitId != null) {
                    System.out.println("No commit with that id exists.");
                    return;
                }
                completeCommitId = completeCommitIdCandidate;
            }
        }
        if (completeCommitId == null) {
            System.out.println("No commit with that id exists.");
            return;
        }
        commitID = completeCommitId;

        if (!getCommitNode(commitID).getBlobs().containsKey(filename)) {
            System.out.println("File does not exist in that commit.");
            return;
        }
        Blob checkfile = getBlob(getCommitNode(commitID).
                getBlobs().get(filename));
        checkfile.checkout();
    }

    public void checkout(String branchName) {
        if (!currentState.containsBranch(branchName)) {
            System.out.println("No such branch exists.");
            return;
        }
        if (currentState.getCurrentBranchTitle().equals(branchName)) {
            System.out.println("No need to checkout the current branch.");
            return;
        }
        currentState.setCurrentBranch(branchName);
        reset(currentState.getBranchHeadCommitNode(branchName));
    }

    public void find(String commitMessage) {
        Boolean error = true;
        ArrayList<CommitNode> commitNodeList = getAllCommitNode();
        for (CommitNode node : commitNodeList) {
            if (node.getMessage().equals(commitMessage)) {
                System.out.println(node.getHashcode());
                error = false;
            }
        }
        if (error) {
            System.out.println("Found no commit with that message.");
        }
    }

    public void status() {
        System.out.println("=== Branches ===");
        TreeSet<String> branchNames = new TreeSet<>(currentState.getAllBranchNames());
        for (String s : branchNames) {
            if (currentState.getCurrentBranchTitle().equals(s)) {
                System.out.print('*');
            }
            System.out.println(s);
        }
        System.out.print("\n");
        System.out.println("=== Staged Files ===");
        stagedset = new TreeSet<>();
        Set<String> blobKey = currentState.getBlobs().keySet();
        for (String s : blobKey) {
            if (!getCurrentCommit().getBlobs().keySet().contains(s)) {
                stagedset.add(s);
            } else if (!currentState.getBlobFilename(s).equals(
                    getCurrentCommit().getBlobFilename(s))) {
                stagedset.add(s);
            }
        }
        for (String s : stagedset) {
            System.out.println(s);
        }
        System.out.print("\n");
        System.out.println("=== Removed Files ===");
        removedset = new TreeSet<>();
        for (String x : getCurrentCommit().getBlobs().keySet()) {
            if (!blobKey.contains(x)) {
                removedset.add(x);
            }
        }
        for (String s : removedset) {
            System.out.println(s);
        }
        System.out.print("\n");
        System.out.println("=== Modifications Not Staged For Commit ===");
        TreeSet<String> filename = new TreeSet<>();
        for (String s : blobKey) {
            if (!new File(s).exists()) {
                filename.add(s + " (deleted)");
            } else if (!Blob.createBlobFromFile(s).
                    getBlobFilename().equals(currentState.getBlobFilename(s))) {
                filename.add(s + " (modified)");
            }
        }
        for (String each : filename) {
            System.out.println(each);
        }
        System.out.print("\n");
        System.out.println("=== Untracked Files ===");
        TreeSet<String> untracked = new TreeSet<>();
        List<String> files = Utils.plainFilenamesIn(".");
        for (String s : files) {
            if (!blobKey.contains(s)
                    && !getCurrentCommit().getBlobs().containsKey(s)) {
                untracked.add(s);
            }
        }
        for (String each : untracked) {
            System.out.println(each);
        }
        System.out.print("\n");
    }


    public static void main(String... args) {
        Main main = new Main();
        if (args.length == 0) {
            return;
        }
        if (args[0].equals("init")) {
            main.init();
        } else {
            if (!new File(".gitlet").exists()) {
                System.out.println("Not in an initialized gitlet directory.");
                return;
            }
            main.loadCurrentState();
            if (args[0].equals("add")) {
                main.add(args[1]);
            } else if (args[0].equals("commit")) {
                main.commit(args[1]);
            } else if (args[0].equals("log")) {
                main.log();
            } else if (args[0].equals("global-log")) {
                main.globalLog();
            } else if (args[0].equals("reset")) {
                main.reset(args[1]);
            } else if (args[0].equals("rm")) {
                main.remove(args[1]);
            } else if (args[0].equals("branch")) {
                main.branch(args[1]);
            } else if (args[0].equals("rm-branch")) {
                main.rmBranch(args[1]);
            } else if (args[0].equals("find")) {
                main.find(args[1]);
            } else if (args[0].equals("status")) {
                main.status();
            } else if (args[0].equals("checkout")) {
                if (args.length == 4) {
                    main.checkout(args[1], args[2], args[3]);
                } else if (args.length == 3) {
                    main.checkout(args[1], args[2]);
                } else {
                    main.checkout(args[1]);
                }
            } else if (args[0].equals("merge")) {
                main.merge(args[1]);
            }
        }
        main.saveCurrentState();
    }
}
