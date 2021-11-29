/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class GFFEntry {

    // constants
    private static final String FEATURE_ID_TAG = "ID=";
    private static final char ATTRIBUTES_SEPERATOR = ';';
    public static final String FASTA_TAG = "fasta_tag=";

    private String sequenceName;
    private String feature;
    private int start;
    private int end;
    private char strand;
    private String attrbutes;
    private boolean modified; // flag to indicate if the GFF entry has been modified since creation
    private boolean valid; // flag to indicate if the GFF entry is valid within the current context

    public GFFEntry(String sequenceName, int start, int end, String attributes) {
        this.sequenceName = sequenceName;
        this.start = start;
        this.end = end;
        this.attrbutes = attributes;
        this.feature = "gene";
        this.strand = '+';
        this.modified = false;
        this.valid = true;
    }

    GFFEntry(String sequenceName, int start, int end, String attributes, String feature, char strand) {
        this(sequenceName, start, end, attributes);
        this.feature = feature;
        this.strand = strand;
    }

    
    /**
     * @return the sequenceName
     */
    public String getSequenceName() {
        return sequenceName;
    }

    /**
     * @param sequenceName the sequenceName to set
     */
    public void setSequenceName(String sequenceName) {
        this.sequenceName = sequenceName;
    }

    /**
     * @return the feature
     */
    public String getFeature() {
        return feature;
    }

    /**
     * @param feature the feature to set
     */
    public void setFeature(String feature) {
        this.feature = feature;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * @return the strand
     */
    public char getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(char strand) {
        this.strand = strand;
    }

    /**
     * @return the attributes
     */
    public String getAttributes() {
        return attrbutes;
    }

    /**
     * @param attributes the attributes to set
     */
    public void setAttributes(String attributes) {
        this.attrbutes = attributes;
    }

    /**
     * @return the modified
     */
    public boolean isModified() {
        return modified;
    }

    /**
     * @param modified the modified to set
     */
    public void setModified(boolean modified) {
        this.modified = modified;
    }

    /**
     * @return the valid
     */
    public boolean isValid() {
        return valid;
    }

    /**
     * @param valid the valid to set
     */
    public void setValid(boolean valid) {
        this.valid = valid;
    }

    void update_start(int start) {
        update_start(start, true);
    }

    public void update_start(int start, boolean set_modified_flag) {
        if (this.start == start) {
            return;
        }
        this.start = start;
        if (set_modified_flag) {
            this.modified = true;
        }
    }

    public void update_end(int end) {
        update_end(end, true);
    }

    public void update_end(int end, boolean set_modified_flag) {
        if (this.end == end) {
            return;
        }
        this.end = end;
        if (set_modified_flag) {
            this.modified = true;
        }
    }

    public String getID() {
        return sequenceName + ":" + feature + ":" + Integer.toString(start) + ":" + Integer.toString(end);
    }

    @Override
    public String toString() {
        return sequenceName + "\t" + getFeatureId() + "\t" + Integer.toString(start) + "\t" + Integer.toString(end);
    }

    public String toStringHeader() {
        // column names should correspond to to_string_short()
        return "Sequence Name\tFeatureID\tStart\tEnd";
    }

    public String getFeatureId() {
        int startPos = attrbutes.indexOf(FEATURE_ID_TAG, 0);
        if (startPos < 0) {
            return "";
        }
        int end_pos = attrbutes.indexOf(ATTRIBUTES_SEPERATOR, startPos);
        if (end_pos < 0) {
            end_pos = attrbutes.length();
        }

        return attrbutes.substring(startPos + FEATURE_ID_TAG.length(), end_pos);

    }

    public static String addIDAttribute(String id, String fastaTag ) {

        StringBuilder sb = new StringBuilder(FEATURE_ID_TAG);
        sb.append(id);
        sb.append(";");
        sb.append(FASTA_TAG);
        sb.append(fastaTag);

        return sb.toString();

    }
    
    public String getEntryString() {
        StringBuilder sb = new StringBuilder(this.sequenceName);
        sb.append("\t.\t");
        sb.append(this.feature);
        sb.append("\t");
        sb.append(Integer.toString(start));
        sb.append("\t");
        sb.append(Integer.toString(end));
        sb.append("\t.\t");
        sb.append(this.strand);
        sb.append("\t.\t");
        sb.append(this.attrbutes);
        
        return sb.toString();

    }

}
