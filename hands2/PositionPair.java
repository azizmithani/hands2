/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hands2;

/**
 *
 * @author Aziz Mithani <aziz.mithani@lums.edu.pk>
 */
public class PositionPair implements Comparable<Object> {

    private int position1;
    private int position2;

    public PositionPair(int position1, int position2) {
        this.position1 = position1;
        this.position2 = position2;
    }

    @Override
    public int compareTo(Object object) {
        PositionPair positionPair = (PositionPair) object;

        if (this.position1 < positionPair.position1) {
            return -1;
        } else if (this.position1 == positionPair.position1) {
            if (this.position2 < positionPair.position2) {
                return -1;
            } else if (this.position2  == positionPair.position2) {
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
        if (!(object instanceof PositionPair)) {
            return false;
        }
        PositionPair positionPair = (PositionPair) object;
        return this.position1 == positionPair.position1 && this.position2 == positionPair.position2;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 31 * hash + this.position1;
        hash = 31 * hash + this.position2;
        return hash;
    }


    @Override
    public String toString() {
        return "(" + position1 + "," + position2 + ")";
    }

    public int getPosition1() {
        return position1;
    }

    public void setPosition1(int position1) {
        this.position1 = position1;
    }

    public int getPosition2() {
        return position2;
    }

    public void setPosition2(int position2) {
        this.position2 = position2;
    }



}