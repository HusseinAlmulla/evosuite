package org.evosuite.ga.metaheuristics;


/*
Original source is in Python
https://github.com/ShangtongZhang/reinforcement-learning-an-introduction
*/

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;

import antlr.collections.List;
import org.evosuite.utils.LoggingUtils;
public class IHT {

	int size;
	int overFullCount;
	public Hashtable  dictionary = new Hashtable();
	public IHT() {
		super();
	}
	
	public IHT(int size) {
		super();
		this.size = size;
		this.overFullCount = 0;
		this.dictionary =  new Hashtable();;
	}

	public int count() {
		return dictionary.size();
	}
	
	public boolean fullp() {
		if (dictionary.size() > size)
			return true;
		else
			return false;
	}
	
	Map<Integer, int[]> index_features = new HashMap<Integer, int[]>();
	
	public int getIndex(int[] obj ) {
		
		Hashtable<int[], Integer> d = dictionary;
        //if (d.contains(obj))
		if(d.get(obj) != null)
        		return (int) d.get(obj);
        
        int sizetmp = size;
        int count = count();
        
        if (count >= sizetmp) {
            if (overFullCount == 0)
            System.out.println("IHT full, starting to allow collisions");
            overFullCount += 1;
            
            String st="";
    		for(int i: obj)
    			st += i;
            int x1 = Math.abs(st.hashCode());
            int x2 = x1 % size; 
            //LoggingUtils.getEvoLogger().info("HASH " + x1 + " size  " + size +  " X2 "  + x2);
            return  x2;
        }
        else {
        	d.put(obj, count);
        	saveFeatureWithIndex(obj);
        	//LoggingUtils.getEvoLogger().info("count " + count);
            return count;
        }
	}
	
	public void saveFeatureWithIndex(int[] obj) {
		Hashtable<int[], Integer> d = dictionary;
		int index = -2;
		if(d.get(obj) != null) {
        	index = (int) d.get(obj);
        	index_features.put(index, obj);
		}
		else
			LoggingUtils.getEvoLogger().info("\nFeature index not found  XXXXXXXXXXXXXXXXXX\n");
	}

	public int hashcoords(ArrayList<Integer> coordinates) {
		int[] arr = coordinates.stream().mapToInt(Integer::intValue).toArray();
	     return getIndex(arr);
	}

	public ArrayList<Integer> tiles(int numTiling, double[] floats, int ints) {
		int[] qfloats = new int[floats.length];

		for(int i=0; i<floats.length; i++) 
			qfloats[i] = (int) Math.floor(floats[i]*numTiling);
		
		ArrayList<Integer> Tiles = new ArrayList<Integer>();
		
		for(int tiling=0; tiling < numTiling; tiling++) {
			int tilingX2 = tiling * 2;
			ArrayList<Integer> coords = new ArrayList<Integer>();
			coords.add(tiling);
			int b = tiling;
			
			for(int q=0; q < qfloats.length; q++) {
				coords.add((int) (Math.floor(q+b)/numTiling));
				b += tilingX2;
			}
			
			coords.add(ints);		
			Tiles.add(hashcoords(coords));
		}
		
		return Tiles;
	}
	
	
	
	
}
