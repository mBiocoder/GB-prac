package com.company;

public class Path {
    private String path;
    private int length;

    public Path() {
        this.path = "";
        this.length = 0;
    }

    public Path(String path, int length) {
        this.path = path;
        this.length = length;
    }

    public void addNodeFront(String node) {
        this.length++;
        if (this.path.length() == 0) this.path = "" + node;
        else this.path = node + "|" + this.path;
    }

    public void addNodeBack(String node) {
        this.length++;
        if (this.path.length() == 0) this.path = "" + node;
        else this.path = this.path + "|" + node;
    }

    public String getPath() {
        return this.path;
    }

    public int getLength() {
        return this.length;
    }

    public Path copy() {
        return new Path(this.path + "", this.length);
    }
}
