/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hands2;

import utility.Base;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class BasePair implements Comparable<Object> {

    private Base firstBase;
    private Base secondBase;
 
    @Override
    public int compareTo(Object object) {
        BasePair snpPair = (BasePair) object;

        if (this.firstBase.compareTo(snpPair.firstBase) < 0) {
            return -1;
        } else if (this.firstBase.compareTo(snpPair.firstBase) == 0) {
            if (this.secondBase.compareTo(snpPair.secondBase) < 0) {
                return -1;
            } else if (this.secondBase.compareTo(snpPair.secondBase) == 0) {
                return 0;
            } else {
                return 1;
            }
        } else {
            return 1;
        }
    }

    @Override
    public boolean equals(Object object) {
        if (!(object instanceof BasePair)) {
            return false;
        }
        BasePair snpPair = (BasePair) object;
        return this.firstBase.equals(snpPair.firstBase) && this.secondBase.equals(snpPair.secondBase);
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 97 * hash + (this.firstBase != null ? this.firstBase.hashCode() : 0);
        hash = 97 * hash + (this.secondBase != null ? this.secondBase.hashCode() : 0);
        return hash;
    }

    public BasePair(Base snp1, Base snp2) {
        this.firstBase = snp1;
        this.secondBase = snp2;
    }

    @Override
    public String toString() {
        return "(" + firstBase + "," + secondBase + ")";
    }

    public Base getFirstBase() {
        return firstBase;
    }

    public Base getSecondBase() {
        return secondBase;
    }


}