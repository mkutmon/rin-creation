package cytargetlinker.conversion.graph;

import java.io.IOException;
import java.io.PrintWriter;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;

/**
 * 
 * @author Thomas
 *
 */
public class XGMMLWriter {
	final static String NS = "http://www.cs.rpi.edu/XGMML";
	
	public static <N, E> void write(Graph graph, PrintWriter out) throws IOException {
		Document doc = new Document();

		Element root = new Element("graph", NS);
		doc.setRootElement(root);
		
		root.setAttribute("id", "" + System.currentTimeMillis());
		root.setAttribute("label", graph.getTitle());
		
		//Print the graph attributes
		printAttributes(graph, root);
		
		//Create the nodes
		for(Node n : graph.getNodes()) {
			Element e = new Element("node");
			e.setAttribute("id", n.getId());
			e.setAttribute("label", n.getId());

			printAttributes(n, e);
			
			root.addContent(e);
		}
		
		//Create the edges
		for(Edge edge : graph.getEdges()) {
			Node src = edge.getSrc();
			Node tgt = edge.getTgt();
			
			Element e = new Element("edge");
			e.setAttribute("id", edge.getId());
			e.setAttribute("label", edge.getId());
			e.setAttribute("source", src.getId());
			e.setAttribute("target", tgt.getId());
			
			Object interaction = edge.getAttribute("interaction");
			Element ie = new Element("att");
			
			ie.setAttribute("label", "interaction");
			ie.setAttribute("name", "interaction");
			ie.setAttribute("value", interaction == null ? "" : interaction.toString());	
			ie.setAttribute("type", "string");
			
			e.addContent(ie);
			
			printAttributes(edge, e);
			
			root.addContent(e);
		}
		
		XMLOutputter xmlcode = new XMLOutputter(Format.getPrettyFormat());
		Format f = xmlcode.getFormat();
		f.setEncoding("UTF-8");
		f.setTextMode(Format.TextMode.PRESERVE);
		xmlcode.setFormat(f);
		xmlcode.output(doc, out);
	}
	
	private static void printAttributes(AttributeHolder attr, Element elm) {
		for(String a : attr.getAttributeNames()) {
			Object o = attr.getAttribute(a);
			if(o == null) continue;
			
			String type =  o instanceof Number ? "real" : "string";
			
			Element e = new Element("att");
			if(o.toString().contains("[") && o.toString().contains("]")) {
				type = "list";
				e.setAttribute("type", "list");
				e.setAttribute("name", "identifiers");
				String [] str = ((String)o).replace("[", "").replace("]", "").split(",");
				for(int i = 0; i < str.length; i++) {
					Element child = new Element("att");
					child.setAttribute("type", "string");
					child.setAttribute("name", "identifiers");
					child.setAttribute("value", str[i]);
					e.addContent(child);
				}
			} else {
				e.setAttribute("label", a);
				e.setAttribute("name", a);
				e.setAttribute("value", "" + o);
				if(a.equals("context+ score") || a.equals("score") || a.equals("pvalue")) {
					e.setAttribute("type", "real");
				} else {
					e.setAttribute("type", type);
				}
			}

			elm.addContent(e);
		}
	}
}