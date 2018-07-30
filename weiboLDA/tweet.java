package weiboLDA;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import Common.Stopwords;

import Common.ComUtil;
import Common.Stopwords;

//import edu.mit.jwi.item.Word;

public class tweet {

	protected int time;
	protected String id;
	protected int[] tweetwords;
	protected int topic;
	public tweet(String dataline,int topic, HashMap<String, Integer> wordMap,
			ArrayList<String> uniWordMap) {

		int number = wordMap.size();
		this.topic=topic;
		// String inline = dataline.substring(20);
		String inline = dataline; // no specific restriction of input data

		ArrayList<Integer> words = new ArrayList<Integer>();
		ArrayList<String> tokens = new ArrayList<String>();
		String []tmp=inline.split("\\s");
		id=tmp[0];
		String []token=tmp[3].split(",");
		
		for (int i = 0; i < token.length; i++) {
			String tmpToken = token[i].trim();
			if (!(tmpToken.isEmpty())&&!Stopwords.isStopword(tmpToken) && !isNoisy(tmpToken)) {
				if (!wordMap.containsKey(tmpToken)) {
					System.out.println(tmpToken);
					words.add(number);
					wordMap.put(tmpToken, number++);
					uniWordMap.add(tmpToken);
				} else {
					words.add(wordMap.get(tmpToken));
				}
			}
		}

		tweetwords = new int[words.size()];

		for (int w = 0; w < words.size(); w++) {
			tweetwords[w] = words.get(w);
		}
		words.clear();

	}

	private boolean isNoisy(String token) {
		if (token.compareTo("​")==0) {
			return true;
		}
		if (token.toLowerCase().contains("#pb#")
				|| token.toLowerCase().contains("http:"))
			return true;
		if (token.contains("@")) {
			return true;
		}
		if (token.contains("#")) {
			return true;
		}
		if (isMessyCode(token)) {
			return true;
		}
		return token.matches("[\\p{Punct}]+");
		// // at least contains a word
		// Pattern MY_PATTERN = Pattern.compile(".*[a-z]+.*");
		// Matcher m = MY_PATTERN.matcher(token.toLowerCase());
		// // filter @xxx
		// if (!m.matches() || token.contains("@")) {
		// return true;
		// }
		// else
		// return false;
	}

    public static boolean isMessyCode(String str) {
        Pattern p = Pattern.compile("\\s*|t*|r*|n*");
        Matcher m = p.matcher(str);
        String after = m.replaceAll("");
        String temp = after.replaceAll("\\p{P}", "");
        char[] ch = temp.trim().toCharArray();
        float chLength = ch.length;
        float count = 0;
        for (int i = 0; i < ch.length; i++) {
            char c = ch[i];
            if (!Character.isLetterOrDigit(c)) {
                if (!isChinese(c)) {
                    count = count + 1;
                }
            }
        }
        float result = count / chLength;
        if (result > 0.4) {
            return true;
        } else {
            return false;
        }
 
    }

    public static boolean isChinese(char c) {
        Character.UnicodeBlock ub = Character.UnicodeBlock.of(c);
        if (ub == Character.UnicodeBlock.CJK_UNIFIED_IDEOGRAPHS
                || ub == Character.UnicodeBlock.CJK_COMPATIBILITY_IDEOGRAPHS
                || ub == Character.UnicodeBlock.CJK_UNIFIED_IDEOGRAPHS_EXTENSION_A
                || ub == Character.UnicodeBlock.GENERAL_PUNCTUATION
                || ub == Character.UnicodeBlock.CJK_SYMBOLS_AND_PUNCTUATION
                || ub == Character.UnicodeBlock.HALFWIDTH_AND_FULLWIDTH_FORMS) {
            return true;
        }
        return false;
    }

}
