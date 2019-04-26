# Module for the Heliophysics Events Knowledgebase
# Contributors: Nicholas Padmanabhan, Mark Cheung

import argparse
import datetime
import numpy as np
import pandas as pd
import xml.etree.cElementTree as ET

import pprint
pp = pprint.PrettyPrinter(indent = 4)

class HER_Event:

	def __init__(self, id):
		self.eventType = id
		self.specFile = "/ssw/vobs/ontology/data/VOEvent_Spec.txt"
		self.Reference_Names = np.zeros(shape = 20, dtype = str)
		self.Reference_Links = np.zeros(shape = 20, dtype = str)
		self.Reference_Types = np.zeros(shape = 20, dtype = str)
		self.Description = ""
		self.Citations = np.zeros(shape = 20, dtype = str)

		self.__voevent_spec = pd.read_csv(self.specFile, skiprows = 2)
		self.__vals = self.__voevent_spec[self.eventType]
		self.__params = self.__voevent_spec["Parameter"].str.upper()
		self.__categories = self.__voevent_spec["VOParamType"]
		self.__sources = self.__voevent_spec["Source"]
		self.__r_o = self.__voevent_spec["R/O"]
		self.__types = self.__voevent_spec["Type"]
		self.__datasetup()

	def __datasetup(self):
		self.required = {}
		self.optional = {}

		for i in range(1, len(self.__vals)):
			if self.__vals[i] == "9":
				if self.__types[i] == "string":
					self.required[self.__params[i]] = "blank"
				elif self.__types[i] == "long":
					self.required[self.__params[i]] = -999999
				elif self.__types[i] == "float":
					self.required[self.__params[i]] = float("inf")
				elif self.__types[i] == "integer":
					self.required[self.__params[i]] = -9999
			elif self.__vals[i] == "5":
				if self.__types[i] == "string":
					self.optional[self.__params[i]] = "blank"
				elif self.__types[i] == "long":
					self.optional[self.__params[i]] = -999999
				elif self.__types[i] == "float":
					self.optional[self.__params[i]] = float("inf")
				elif self.__types[i] == "integer":
					self.optional[self.__params[i]] = -9999

		defaulted_req_keys = ["EVENT_TYPE", "KB_ARCHIVDATE", "KB_ARCHIVID", 
							  "KB_ARCHIVIST", "KB_ARCHIVURL", "EVENT_COORDSYS", 
							  "EVENT_ENDTIME", "EVENT_STARTTIME"]

		defaulted_req_values = ["{}: {}".format(self.eventType, self.__vals[0]),
								"Reserved for KB archivist: KB entry date",
								"Reserved for KB archivist: KB entry identifier",
								"Reserved for KB archivist: KB entry made by",
								"Reserved for KB archivist: URL to suppl. info.",
								"UTC-HPC-TOPO", "1492-10-12 00:00:00", "1492-10-12 00:00:00"]

		for i in range(len(defaulted_req_keys)):
			self.required[defaulted_req_keys[i]] = defaulted_req_values[i]

		defaulted_opt_keys = ["EVENT_EXPIRES"]

		defaulted_opt_values = ["1492-10-12 00:00:00"]

		for i in range(len(defaulted_opt_keys)):
			self.optional[defaulted_opt_keys[i]] = defaulted_opt_values[i]

	def exportEvent(self, filename = None, suffix = ""):
		remove_chars = [' ',':',';','/','\\', '{{','}}' ,'[[',']]', '~', '=', '.', ',', '-', '+']
		for r in remove_chars:
			self.required["FRM_NAME"] = self.required["FRM_NAME"].replace(r,"")

		d = datetime.datetime.utcnow().isoformat()
		eventIdentifier = (self.eventType.split(':'))[0] + "_" + self.required["FRM_NAME"] + "_" + d + suffix
		self.required["KB_ARCHIVID"] = "ivo://helio-informatics.org/{}".format(eventIdentifier)
		print(self.required["KB_ARCHIVID"])

		""" Sets up the XML data tree using Python's cElementTree library """

		voe = ET.Element("voe:VOEvent")
		voe_info = ET.Comment("This is a VOEvent, generated from the annotation process")
		voe.insert(1, voe_info)
		voe.set("xmlns:voe", "http://www.ivoa.net/xml/VOEvent/v1.1")
		voe.set("xmlns:stc", "http://www.ivoa.net/xml/STC/stc-v1.30.xsd")
		voe.set("xmlns:lmsal", "http://www.lmsal.com/helio-informatics/lmsal-v1.0.xsd")
		voe.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
		voe.set("xmlns:crd", "urn:nvo-coords")
		voe.set("role", "observation")
		voe.set("version", "1.1")
		voe.set("xsi:schemaLocation", "http://www.ivoa.net/xml/VOEvent/v1.1 http://www.lmsal.com/helio-informatics/VOEvent-v1.1.xsd")
		voe.set("ivorn", self.required["KB_ARCHIVID"])

		""" Who section """

		Who = ET.SubElement(voe, "Who")
		who_info = ET.Comment("Data pertaining to curation")
		Who.insert(1, who_info)

		AuthorIVORN = ET.SubElement(Who, "AuthorIVORN")
		Author = ET.SubElement(Who, "Author")
		contactName = ET.SubElement(Author, "contactName")
		contactEmail = ET.SubElement(Author, "contactEmail")
		Date = ET.SubElement(Who, "Date")

		AuthorIVORN.text = self.required["KB_ARCHIVURL"]
		contactName.text = self.required["FRM_NAME"]
		contactEmail.text = self.required["FRM_CONTACT"]
		Date.text = d

		""" What section """

		What = ET.SubElement(voe, "What")
		what_info = ET.Comment("Data about what was measured/observed")
		What.insert(1, what_info)

		Description = ET.SubElement(What, "Description")
		Description.text = self.Description

		What_Group_Required = ET.SubElement(What, "Group")
		What_Group_Required.set("name", "{}_required".format(self.__vals[0]))

		What_Group_Optional = ET.SubElement(What, "Group")
		What_Group_Optional.set("name", "{}_optional".format(self.__vals[0]))

		for param in self.required:
			i = pd.Index(self.__params).get_loc(param)
			if self.__categories[i] == "what":
				t = ET.SubElement(What_Group_Required, "Param")
				t.set("name", param)
				t.set("value", str(self.required[param]))

		for param in self.optional:
			i = pd.Index(self.__params).get_loc(param)
			if self.optional[param] != float("inf") and self.optional[param] != "blank" and self.optional[param] != -9999 and self.optional[param] != -999999 and self.optional[param] != "1492-10-12 00:00:00":
				if self.__categories[i] == "what":
					t = ET.SubElement(What_Group_Optional, "Param")
					t.set("name", param)
					t.set("value", str(self.optional[param]))

		""" Where-When section """

		WhereWhen = ET.SubElement(voe, "WhereWhen")
		where_when_info = ET.Comment("Data pertaining to where and when something occurred")
		WhereWhen.insert(1, where_when_info)

		ObsDataLocation = ET.SubElement(WhereWhen, "ObsDataLocation")
		ObsDataLocation.set("xmlns", "http://www.ivoa.net/xml/STC/stc-v1.30.xsd")

		ObservatoryLocation = ET.SubElement(ObsDataLocation, "ObservatoryLocation")
		ObservatoryLocation_AstroCoordSystem = ET.SubElement(ObservatoryLocation, "AstroCoordSystem")

		ObservatoryLocation_AstroCoords = ET.SubElement(ObservatoryLocation, "AstroCoords")

		ObservationLocation = ET.SubElement(ObsDataLocation, "ObservationLocation")
		ObservationLocation_AstroCoordSystem = ET.SubElement(ObservationLocation, "AstroCoordSystem")
		ObservationLocation_AstroCoords = ET.SubElement(ObservationLocation, "AstroCoords")

		Time = ET.SubElement(ObservationLocation_AstroCoords, "Time")
		TimeInstant = ET.SubElement(Time, "TimeInstant")
		TimeInstant_ISOTime = ET.SubElement(TimeInstant, "ISOTime")
		Position2D = ET.SubElement(ObservationLocation_AstroCoords, "Position2D")
		Value2 = ET.SubElement(Position2D, "Value2")
		Position2D_Value_C1 = ET.SubElement(Value2, "C1")
		Position2D_Value_C2 = ET.SubElement(Value2, "C2")
		Error2 = ET.SubElement(Position2D, "Error2")
		Position2D_Error_C1 = ET.SubElement(Error2, "C1")
		Position2D_Error_C2 = ET.SubElement(Error2, "C2")
		AstroCoordArea = ET.SubElement(ObservationLocation, "AstroCoordArea")
		TimeInterval = ET.SubElement(AstroCoordArea, "TimeInterval")
		StartTime = ET.SubElement(TimeInterval, "StartTime")
		StartTime_ISOTime = ET.SubElement(StartTime, "ISOTime")
		StopTime = ET.SubElement(TimeInterval, "StopTime")
		StopTime_ISOTime = ET.SubElement(StopTime, "ISOTime")
		Box2 = ET.SubElement(AstroCoordArea, "Box2")
		Center = ET.SubElement(Box2, "Center")
		Box_Center_C1 = ET.SubElement(Center, "C1")
		Box_Center_C2 = ET.SubElement(Center, "C2")
		Size = ET.SubElement(Box2, "Size")
		Box_Size_C1 = ET.SubElement(Size, "C1")
		Box_Size_C2 = ET.SubElement(Size, "C2")

		WhereWhen_Group_Required = ET.SubElement(WhereWhen, "Group")
		WhereWhen_Group_Required.set("name", "{}_required".format(self.__vals[0]))

		WhereWhen_Group_Optional = ET.SubElement(WhereWhen, "Group")
		WhereWhen_Group_Optional.set("name", "{}_optional".format(self.__vals[0]))

		for param in self.required:
			i = pd.Index(self.__params).get_loc(param)
			if self.__categories[i] == "wherewhen" and self.__sources[i] == "data": ############# HERE DATA?? EVENT PEAKTIME?? ########
				t = ET.SubElement(WhereWhen_Group_Required, "Param")
				t.set("name", param)
				t.set("value", str(self.required[param]))

		for param in self.optional:
			i = pd.Index(self.__params).get_loc(param)
			if self.optional[param] != float("inf") and self.optional[param] != "blank" and self.optional[param] != -9999 and self.optional[param] != -999999 and self.optional[param] != "1492-10-12 00:00:00":
				if self.__categories[i] == "wherewhen":
					t = ET.SubElement(WhereWhen_Group_Optional, "Param")
					t.set("name", param)
					t.set("value", str(self.optional[param]))

		ObservationLocation.set("id", self.required["OBS_OBSERVATORY"])
		TimeInstant_ISOTime.text = str(self.required["EVENT_STARTTIME"])
		Position2D.set("unit", self.required["EVENT_COORDUNIT"])
		Position2D_Value_C1.text = str(self.required["EVENT_COORD1"])
		Position2D_Value_C2.text = str(self.required["EVENT_COORD2"])
		Position2D_Error_C1.text = str(self.required["EVENT_C1ERROR"])
		Position2D_Error_C2.text = str(self.required["EVENT_C2ERROR"])
		StartTime_ISOTime.text = str(self.required["EVENT_STARTTIME"])
		StopTime_ISOTime.text = str(self.required["EVENT_ENDTIME"])

		"""
		The user-inputted data requires the lower-left and upper-right coordinates of the bounding box.
		The XML output contains the bounding box center and size.
		The below four lines use the lower-left and upper-right coordinates to calculate the center and size of the bounding box.
		"""
		Box_Center_C1.text = "%f" % ((float(self.required["BOUNDBOX_C1LL"]) + float(self.required["BOUNDBOX_C1UR"])) / 2.0)
		Box_Center_C2.text = "%f" % ((float(self.required["BOUNDBOX_C2LL"]) + float(self.required["BOUNDBOX_C2UR"])) / 2.0)
		Box_Size_C1.text = "%f" % np.abs((float(self.required["BOUNDBOX_C1LL"]) - float(self.required["BOUNDBOX_C1UR"])))
		Box_Size_C2.text = "%f" % np.abs((float(self.required["BOUNDBOX_C2LL"]) - float(self.required["BOUNDBOX_C2UR"])))

		ObservatoryLocation_AstroCoords.set("id", self.required["EVENT_COORDSYS"])
		ObservatoryLocation_AstroCoords.set("coord_system_id", self.required["EVENT_COORDSYS"])
		ObservationLocation_AstroCoords.set("coord_system_id", self.required["EVENT_COORDSYS"])

		""" How section """

		How = ET.SubElement(voe, "How")
		how_info = ET.Comment("Data pertaining to how the feature/event detection was performed")
		How.insert(1, how_info)

		lmsal_data = ET.SubElement(How, "lmsal:data")
		lmsal_method = ET.SubElement(How, "lmsal:method")

		How_Group_Optional = ET.SubElement(How, "Group")
		How_Group_Optional.set("name", "{}_optional".format(self.__vals[0]))

		for param in self.required:
			i = pd.Index(self.__params).get_loc(param)
			if self.__categories[i] == "how" and self.__sources[i] == "method":
				if param != "FRM_URL":
					t = ET.SubElement(lmsal_method, "lmsal:{}".format(param))
					t.text = str(self.required[param])

			if self.__categories[i] == "how" and self.__sources[i] == "data":
				t = ET.SubElement(lmsal_data, "lmsal:{}".format(param))
				t.text = str(self.required[param])

		for param in self.optional:
			i = pd.Index(self.__params).get_loc(param)
			if self.optional[param] != float("inf") and self.optional[param] != "blank" and self.optional[param] != -9999 and self.optional[param] != -999999 and self.optional[param] != "1492-10-12 00:00:00":
				if self.__categories[i] == "how":
					t = ET.SubElement(How_Group_Optional, "Param")
					t.set("name", param)
					t.set("value", str(self.optional[param]))

		""" Why section """

		Why = ET.SubElement(voe, "Why")
		Inference = ET.SubElement(Why, "Inference")
		Concept = ET.SubElement(Why, "Concept")
		lmsal_EVENT_TYPE = ET.SubElement(Why, "lmsal:EVENT_TYPE")

		if self.optional["EVENT_PROBABILITY"] != float("inf"):
			Inference.set("probability", str(self.optional["EVENT_PROBABILITY"]))

		Concept.text = self.__vals[0]
		lmsal_EVENT_TYPE.text = "{}: {}".format(self.eventType, self.__vals[0])

		Why_Group_Required = ET.SubElement(Why, "Group")
		Why_Group_Required.set("name", "{}_required".format(self.__vals[0]))

		Why_Group_Optional = ET.SubElement(Why, "Group")
		Why_Group_Optional.set("name", "{}_optional".format(self.__vals[0]))

		for param in self.required:
			i = pd.Index(self.__params).get_loc(param)
			if self.__categories[i] == "why":
				t = ET.SubElement(Why_Group_Required, "Param")
				t.set("name", param)
				t.set("value", str(self.required[param]))

		for param in self.optional:
			i = pd.Index(self.__params).get_loc(param)
			if self.optional[param] != float("inf") and self.optional[param] != "blank" and self.optional[param] != -9999 and self.optional[param] != -999999 and self.optional[param] != "1492-10-12 00:00:00":
				if self.__categories[i] == "why":
					t = ET.SubElement(Why_Group_Optional, "Param")
					t.set("name", param)
					t.set("value", str(self.optional[param]))

		""" Citations section - omitted for now """

		# Citations = ET.SubElement(voe, "Citations")

		# for citation in self.Citations:
		# 	if citation != "":
		# 		pass

		""" References section """

		for i in range(len(self.Reference_Names)):
			if self.Reference_Names[i] != "":
				Reference = ET.SubElement(voe, "Reference")
				Reference.set("name", self.Reference_Names[i])
				Reference.set("type", self.Reference_Types[i])
				Reference.set("uri", self.Reference_Links[i])

		if self.required["FRM_URL"] == "blank":
			self.required["FRM_URL"] = "n/a"

		Reference = ET.SubElement(voe, "Reference")
		Reference.set("name", "FRM_URL")
		Reference.set("uri", self.required["FRM_URL"])

		if self.optional["OBS_DATAPREPURL"] != "blank":
			Reference = ET.SubElement(voe, "Reference")
			Reference.set("name", "OBS_DATAPREPURL")
			Reference.set("uri", self.optional["OBS_DATAPREPURL"])

		if self.optional["EVENT_MAPURL"] != "blank":
			Reference = ET.SubElement(voe, "Reference")
			Reference.set("name", "EVENT_MAPURL")
			Reference.set("uri", self.optional["EVENT_MAPURL"])

		""" Clean up and export the XML file """

		self.__indentxml(voe)
		tree = ET.ElementTree(voe)

		""" Checks if all required values have been inputted """
		to_add = []
		for param in self.required:
			i = pd.Index(self.__params).get_loc(param)
			if self.required[param] == float("inf") or self.required[param] == "blank" or self.required[param] == -9999 or self.required[param] == -999999 or self.required[param] == "1492-10-12 00:00:00":
				to_add.append(param)
		
		if len(to_add) != 0:
			print("\n*** Error - required value(s) not entered:")
			pp.pprint(to_add)
			return

		""" Removes empty required/optional elements """
		for i in range(len(voe)):
			for group in voe[i].findall("Group"):
				if group.getchildren() == []:
					voe[i].remove(group)

		""" Removes Inference element if no probability is entered """
		if voe[5][0].attrib == {}:
			voe[5].remove(voe[5][0])

		if filename is None:
			filename = eventIdentifier + ".xml"
			print(filename)

		tree.write(filename)

	def __str__(self):
		self.data = {
			"REQUIRED" : self.required,
			"OPTIONAL" : self.optional,
			"SPECFILE" : self.specFile,
			"REFERENCE_NAMES" : self.Reference_Names,
			"REFERENCE_LINKS" : self.Reference_Links,
			"REFERENCE_TYPES" : self.Reference_Types,
			"DESCRIPTION" : self.Description,
			"CITATIONS" : self.Citations
		}
		pp.pprint(self.data)
		return ""

	def __indentxml(self, elem, level = 0):
		i = "\n" + level * "	"
		if len(elem):
			if not elem.text or not elem.text.strip():
				elem.text = i + "	"
			if not elem.tail or not elem.tail.strip():
				elem.tail = i
			for elem in elem:
				self.__indentxml(elem, level + 1)
			if not elem.tail or not elem.tail.strip():
				elem.tail = i
		else:
			if level and (not elem.tail or not elem.tail.strip()):
				elem.tail = i
