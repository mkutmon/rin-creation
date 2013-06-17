package cytargetlinker.conversion.data;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;

public class MTI {

	private MiRNANode source;
	private GeneNode target;
	private String score;
	private String interactionClass = "microRNA-target";
	
	public MTI(MiRNANode source, GeneNode target, String score) {
		super();
		this.source = source;
		this.target = target;
		this.score = score;
	}
	
	public Edge createEdge(Graph graph, String datasource, int count) {
		Edge e = graph.addEdge("" + count, source.getNode(graph), target.getNode(graph));
		e.setAttribute("score", getScore().toString());
		e.setAttribute("interactionType", getInteractionClass());
		e.setAttribute("datasource", datasource);
		return e;
	}
	
	public MiRNANode getSource() {
		return source;
	}
	public void setSource(MiRNANode source) {
		this.source = source;
	}
	public GeneNode getTarget() {
		return target;
	}
	public void setTarget(GeneNode target) {
		this.target = target;
	}
	public String getScore() {
		return score;
	}
	public void setScore(String score) {
		this.score = score;
	}
	public String getInteractionClass() {
		return interactionClass;
	}
	public void setInteractionClass(String interactionClass) {
		this.interactionClass = interactionClass;
	}
}
