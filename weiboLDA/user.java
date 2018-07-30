package weiboLDA;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Common.FileUtil;

//import edu.mit.jwi.item.Word;

public class user {

	protected String userID;
	
	protected int tweetCnt;
	
	public ArrayList<tweet> tweets = new ArrayList<tweet>();
	
	
	public user(String Dir, String id, HashMap<String, Integer> wordMap, 
			ArrayList<String> uniWordMap) {
		
		this.userID = id;
		ArrayList<String> datalines = new ArrayList<String>();
		ArrayList<String> filelist = new ArrayList<String>();
		FileUtil.readLines(Dir, filelist);
		for(int fileNo=0;fileNo<filelist.size();fileNo++) {
			String datadir=filelist.get(fileNo);
		FileUtil.readLines(datadir, datalines);		
		
		this.tweetCnt += datalines.size();
		
		for(int lineNo = 0; lineNo < datalines.size(); lineNo++) {
			String line = datalines.get(lineNo);
			tweet tw = new tweet(line,fileNo, wordMap, uniWordMap);
			tweets.add(tw);						
		}
		datalines.clear();
		}
	}
}
