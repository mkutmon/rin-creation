package cytargetlinker.conversion.graph;

import java.io.PrintWriter;

import org.freehep.util.io.IndentPrintWriter;

import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;

/**
 * 
 * @author Thomas
 *
 */
public class GmlWriter {
	public static void write(Graph graph, PrintWriter out) {
		IndentPrintWriter iout = new IndentPrintWriter(out);
		iout.setIndentString("\t");
		int indent = 0;
		
		iout.println("graph [");
		iout.setIndent(++indent);

		//Print nodes and attributes
		for(Node n : graph.getNodes()) {
			iout.println("node [");
			iout.setIndent(++indent);
			iout.println("id\t" + n.hashCode());
			iout.println("identifier\t" + '"' + n.getId() + '"');
			printAttributes(iout, n);
			iout.setIndent(--indent);
			iout.println("]");
		}
		//Print edges and attributes
		for(Edge e : graph.getEdges()) {
			iout.println("edge [");
			iout.setIndent(++indent);
			String srcS = e.getSrc().getId();
			String tgtS = e.getTgt().getId();
			int src = srcS.hashCode();
			int tgt = tgtS.hashCode();
			iout.println("source\t"  + src);
			iout.println("target\t"  + tgt);
			iout.println("id\t" + '"' + (src > tgt ? tgt + "," + src : src + "," + tgt) + '"');
			iout.println("identifier\t" + '"' + (src > tgt ? tgtS + "," + srcS : srcS + "," + tgtS) + '"');
			printAttributes(iout, e);
			iout.setIndent(--indent);
			iout.println("]");
		}
		//Print network attributes
		printAttributes(iout, graph);
		iout.setIndent(--indent);
		iout.println("]");
	}
	
	private static void printAttributes(PrintWriter out, AttributeHolder attributes) {
		for(String k : attributes.getAttributeNames()) {
			//Find out if this is a number
			Object v = attributes.getAttribute(k);
			if(v == null || "".equals(v)) continue; //Skip empty attributes
			boolean isNumber = v instanceof Number;
			if(!isNumber) {
				String s = v.toString();
				s = s.replaceAll("\"", "");
				v = '"' + s + '"';
			}
			out.println(k + "\t" + v);
		}
	}
}
