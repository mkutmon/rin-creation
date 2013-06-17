package cytargetlinker.conversion.utils;

/**
 * common attributes of a xgmml necessary
 * for the CyTargetLinker plugin for Cytoscape
 * @author Thomas
 *
 */
public enum CommonAttributes {
	DATABASE("Database"),
	DB_Version("DatabaseVersion"),
	SOURCE_FILE("SourceFile"),
//	SOURCE_ID_ATTRIBUTE("SourceIdAttribute"),
//	TARGET_ID_ATTRIBUTE("TargetIdAttribute"),
	SOURCE_DATASOURCE("SourceDatasource"),
	TARGET_DATASOURCE("TargetDatasource"),
	SOURCE_TYPE("SourceType"),
	TARGET_TYPE("TargetType"),
	TYPE("Type"),
	IDENTIFIERS("Identifiers");
	
	String name;
	
	private CommonAttributes(String name) {
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
}
