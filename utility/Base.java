/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

/**
 *
 * @author aziz
 */
public class Base implements Comparable<Object> {
    
    private int position;
    private char nucleotide;

        public Base() {
        position = -1;
        nucleotide = DNA.BASE_BLANK;
    }

    public Base(int position, char base) {
        this.position = position;
        this.nucleotide = base;
  }

    public Base(int position, char base, char refBase) {
        this.position = position;
        this.nucleotide = base;
    }

    public char getNucleotide() {
        return nucleotide;
    }

    public void setNucleotide(char nucleotide) {
        this.nucleotide = nucleotide;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }


    @Override
    public int compareTo(Object object) {
        Base base = (Base) object;
        if (this.position < base.position) {
            return -1;
        } else if (this.position > base.position) {
            return 1;
        } else if (this.nucleotide < base.nucleotide) {
            return -1;
        } else if (this.nucleotide > base.nucleotide) {
            return 1;
        } else {
            return 0;
        }

    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 59 * hash + this.position;
        hash = 59 * hash + this.nucleotide;
        return hash;
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof Base)) {
            return false;
        }
        Base base = (Base) object;
        return this.position == base.position && this.nucleotide == base.nucleotide;
    }

    @Override
    public String toString() {
        return position + ":" + nucleotide;
    }

}
