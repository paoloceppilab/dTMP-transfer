import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const DEFAULT_PREFIX = "kras_g12_g13_hotspot_frequency_tumor_only_nonhdr";
const CATEGORIES = [
  "WT",
  "G12D",
  "G12V",
  "G12C",
  "G12A",
  "G12R",
  "G12S",
  "G13D",
  "G13C",
  "G13R",
  "other",
];

function parseArgs(argv) {
  const args = {
    packageDir: null,
    prefix: DEFAULT_PREFIX,
  };
  for (let i = 0; i < argv.length; i += 1) {
    const arg = argv[i];
    if (arg === "--package-dir") {
      args.packageDir = argv[++i];
    } else if (arg === "--prefix") {
      args.prefix = argv[++i];
    } else {
      throw new Error(`Unknown argument: ${arg}`);
    }
  }
  if (!args.packageDir) {
    throw new Error("Missing required argument: --package-dir");
  }
  return args;
}

function colName(idx0) {
  let n = idx0 + 1;
  let out = "";
  while (n > 0) {
    const rem = (n - 1) % 26;
    out = String.fromCharCode(65 + rem) + out;
    n = Math.floor((n - 1) / 26);
  }
  return out;
}

function rangeA1(row1, col1, rows, cols) {
  const start = `${colName(col1 - 1)}${row1}`;
  const end = `${colName(col1 + cols - 2)}${row1 + rows - 1}`;
  return `${start}:${end}`;
}

function parseValue(value) {
  if (value === "") return "";
  if (value === "True") return true;
  if (value === "False") return false;
  if (/^-?(?:\d+|\d*\.\d+)(?:[eE][+-]?\d+)?$/.test(value)) return Number(value);
  return value;
}

async function readTsv(filePath) {
  const text = await fs.readFile(filePath, "utf8");
  const lines = text.trimEnd().split(/\r?\n/);
  const headers = lines[0].split("\t");
  const rows = lines.slice(1).map((line) => {
    const values = line.split("\t").map(parseValue);
    return Object.fromEntries(headers.map((header, idx) => [header, values[idx] ?? ""]));
  });
  return { headers, rows };
}

async function sha256(filePath) {
  const data = await fs.readFile(filePath);
  return crypto.createHash("sha256").update(data).digest("hex");
}

function pct(count, denominator) {
  return denominator ? (100 * count) / denominator : 0;
}

function rowsToMatrix(headers, rows) {
  return [headers, ...rows.map((row) => headers.map((header) => row[header] ?? ""))];
}

function styleTable(sheet, rows, cols, options = {}) {
  const used = sheet.getRange(rangeA1(1, 1, rows, cols));
  used.format.font = { name: "Calibri", size: 10 };
  used.format.borders = { preset: "outside", style: "thin", color: "#D1D5DB" };
  used.format.wrapText = true;
  const header = sheet.getRange(rangeA1(1, 1, 1, cols));
  header.format.fill = options.headerFill ?? "#1F4E79";
  header.format.font = { color: "#FFFFFF", bold: true };
  header.format.horizontalAlignment = "center";
  sheet.freezePanes.freezeRows(1);
  used.format.autofitColumns();
  used.format.autofitRows();
}

function writeSheet(workbook, name, headers, rows, options = {}) {
  const sheet = workbook.worksheets.add(name);
  const matrix = rowsToMatrix(headers, rows);
  sheet.getRange(rangeA1(1, 1, matrix.length, headers.length)).values = matrix;
  styleTable(sheet, matrix.length, headers.length, options);
  return sheet;
}

function buildReadmeRows(runInfo, workbookPath) {
  return [
    { field: "package_name", value: path.basename(path.dirname(path.dirname(workbookPath))) },
    { field: "created_utc", value: runInfo.created_utc },
    { field: "purpose", value: "Corrected Supplementary Figure 5D KRAS codon 12/13 point-mutation spectrum." },
    { field: "figure_title", value: runInfo.figure_title },
    { field: "denominator", value: runInfo.denominator_type },
    { field: "tumor_filter", value: runInfo.filtering_rule.include_tumor_samples_only },
    { field: "hdr_filter", value: runInfo.filtering_rule.exclude_hdr_barcode_positive },
    { field: "interpretation", value: runInfo.corrected_nonhdr.interpretation_note },
    { field: "source_allele_table", value: runInfo.inputs.alleles },
    { field: "figure_ready_tsv", value: runInfo.outputs.figure_ready_tsv },
    { field: "run_info_json", value: runInfo.outputs.run_info },
  ];
}

function buildAggregateRows(runInfo) {
  const denominator = runInfo.corrected_nonhdr.aggregate_denominator_allele_reads;
  return CATEGORIES.map((category) => ({
    plot_name: "kras_g12_g13_hotspot_frequency_tumor_only_nonhdr",
    denominator_type: runInfo.denominator_type,
    aggregate_denominator_reads: denominator,
    category,
    aggregate_reads: runInfo.corrected_nonhdr.aggregate_counts[category],
    aggregate_percent: runInfo.corrected_nonhdr.aggregate_percentages[category],
  }));
}

function buildLongRows(tsvRows, figureReadyPath, runInfo) {
  return tsvRows.flatMap((row) =>
    CATEGORIES.map((category) => ({
      plot_name: "kras_g12_g13_hotspot_frequency_tumor_only_nonhdr",
      specimen_id: row.specimen_id,
      cohort: row.cohort,
      denominator_type: runInfo.denominator_type,
      full_tumor_allele_reads_before_hdr_filter: row.full_tumor_allele_reads_before_hdr_filter,
      strict_hdr_barcode_yes_reads_excluded: row.strict_hdr_barcode_yes_reads_excluded,
      denominator_reads: row.denominator_nonhdr_tumor_allele_reads,
      category,
      reads: row[`${category}_reads`],
      percent: row[`${category}_pct`],
      category_pct_sum: row.category_pct_sum,
      qc_pass: row.category_count_check_pass,
      source_figure_ready_tsv: figureReadyPath,
    })),
  );
}

function buildRemovedHdrRows(runInfo) {
  const oldDenominator = runInfo.old_bundled_reference.tumor_denominator_allele_reads;
  const correctedDenominator = runInfo.corrected_nonhdr.aggregate_denominator_allele_reads;
  const removedDenominator = runInfo.hdr_positive_removed.tumor_strict_hdr_barcode_yes_reads;
  return CATEGORIES.map((category) => {
    const oldReads = runInfo.old_bundled_reference.aggregate_counts[category] ?? 0;
    const removedReads = runInfo.hdr_positive_removed.aggregate_counts[category] ?? 0;
    const correctedReads = runInfo.corrected_nonhdr.aggregate_counts[category] ?? 0;
    return {
      category,
      old_tumor_only_reads: oldReads,
      old_tumor_only_percent: pct(oldReads, oldDenominator),
      hdr_positive_reads_removed: removedReads,
      hdr_positive_percent_of_old_denominator: pct(removedReads, oldDenominator),
      corrected_nonhdr_reads: correctedReads,
      corrected_nonhdr_percent: pct(correctedReads, correctedDenominator),
      denominator_old_tumor_only_reads: oldDenominator,
      denominator_hdr_positive_removed_reads: removedDenominator,
      denominator_corrected_nonhdr_reads: correctedDenominator,
    };
  });
}

function buildQcRows(runInfo) {
  return Object.entries(runInfo.qc_checks).map(([check, pass]) => ({
    qc_check: check,
    pass,
  }));
}

function buildProvenanceRows(runInfo) {
  const rows = [
    { item_type: "input", item: "alleles", path: runInfo.inputs.alleles, sha256: runInfo.checksums_sha256.inputs.alleles },
  ];
  for (const [name, outputPath] of Object.entries(runInfo.outputs)) {
    if (name === "run_info") continue;
    rows.push({
      item_type: "output",
      item: name,
      path: outputPath,
      sha256: runInfo.checksums_sha256.outputs?.[name] ?? "",
    });
  }
  rows.push({
    item_type: "script",
    item: "plot_script",
    path: runInfo.script,
    sha256: runInfo.checksums_sha256.script,
  });
  if (runInfo.workbook_builder?.script) {
    rows.push({
      item_type: "script",
      item: "workbook_builder",
      path: runInfo.workbook_builder.script,
      sha256: runInfo.checksums_sha256.workbook_builder_script ?? "",
    });
  }
  return rows;
}

function setUsefulWidths(workbook) {
  const widthMap = {
    README: { A: 230, B: 720 },
    aggregate_summary: { A: 320, B: 320, C: 180, D: 90, E: 140, F: 140 },
    g12_g13_nonhdr_data: { A: 320, B: 110, C: 80, D: 320, E: 180, F: 170, G: 150, H: 90, I: 120, J: 110 },
    per_sample_summary: { A: 110, B: 80, C: 80, D: 320, E: 170, F: 150, G: 170 },
    removed_hdr_summary: { A: 90, B: 150, C: 150, D: 150, E: 170, F: 150, G: 150 },
    qc_checks: { A: 420, B: 80 },
    provenance: { A: 100, B: 180, C: 720, D: 520 },
  };
  for (const [sheetName, widths] of Object.entries(widthMap)) {
    const sheet = workbook.worksheets.getItem(sheetName);
    for (const [col, widthPx] of Object.entries(widths)) {
      sheet.getRange(`${col}:${col}`).format.columnWidthPx = widthPx;
    }
  }
}

async function compactVerify(workbook) {
  const summary = await workbook.inspect({
    kind: "table",
    range: "aggregate_summary!A1:F12",
    include: "values,formulas",
    tableMaxRows: 12,
    tableMaxCols: 6,
  });
  console.log(summary.ndjson);
  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 300 },
    summary: "final formula error scan",
  });
  console.log(errors.ndjson);
  for (const sheetName of [
    "README",
    "aggregate_summary",
    "g12_g13_nonhdr_data",
    "per_sample_summary",
    "removed_hdr_summary",
    "qc_checks",
    "provenance",
  ]) {
    await workbook.render({ sheetName, range: "A1:K25", scale: 1 });
  }
}

async function writeManifest(packageDir, runInfo, workbookPath) {
  const manifestPath = path.join(packageDir, "metadata", "package_manifest.tsv");
  const outputEntries = [
    ["figure_png", runInfo.outputs.figure_png],
    ["figure_pdf", runInfo.outputs.figure_pdf],
    ["figure_svg", runInfo.outputs.figure_svg],
    ["figure_ready_tsv", runInfo.outputs.figure_ready_tsv],
    ["workbook_xlsx", workbookPath],
    ["run_info_json", runInfo.outputs.run_info],
    ["readme", runInfo.outputs.readme],
  ];
  const rows = [];
  for (const [artifactKey, filePath] of outputEntries) {
    const stat = await fs.stat(filePath);
    rows.push({
      artifact_type: "corrected_nonhdr_sf5d_output",
      artifact_key: artifactKey,
      relative_path: path.relative(packageDir, filePath),
      path: filePath,
      size_bytes: stat.size,
      sha256: await sha256(filePath),
    });
  }
  const headers = ["artifact_type", "artifact_key", "relative_path", "path", "size_bytes", "sha256"];
  const text = [
    headers.join("\t"),
    ...rows.map((row) => headers.map((header) => row[header]).join("\t")),
  ].join("\n") + "\n";
  await fs.writeFile(manifestPath, text, "utf8");
  return manifestPath;
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const packageDir = path.resolve(args.packageDir);
  const prefix = args.prefix;
  const figureReadyPath = path.join(packageDir, "data", `${prefix}_figure_ready.tsv`);
  const runInfoPath = path.join(packageDir, "metadata", `${prefix}_run_info.json`);
  const workbookPath = path.join(packageDir, "data", `${prefix}_data.xlsx`);
  const builderScript = fileURLToPath(import.meta.url);

  const { rows: tsvRows } = await readTsv(figureReadyPath);
  const runInfo = JSON.parse(await fs.readFile(runInfoPath, "utf8"));
  runInfo.workbook_builder = {
    script: builderScript,
    argv: process.argv,
    node_version: process.version,
    artifact_tool: "SpreadsheetFile.exportXlsx from @oai/artifact-tool",
  };
  runInfo.checksums_sha256.workbook_builder_script = await sha256(builderScript);

  const workbook = Workbook.create();
  let readmeSheet;
  try {
    readmeSheet = workbook.worksheets.getActiveWorksheet();
    readmeSheet.name = "README";
    readmeSheet.reset();
  } catch {
    readmeSheet = workbook.worksheets.add("README");
  }

  const readmeRows = buildReadmeRows(runInfo, workbookPath);
  const readmeMatrix = rowsToMatrix(["field", "value"], readmeRows);
  readmeSheet.getRange(rangeA1(1, 1, readmeMatrix.length, 2)).values = readmeMatrix;
  styleTable(readmeSheet, readmeMatrix.length, 2, { headerFill: "#1F4E79" });

  writeSheet(workbook, "aggregate_summary", Object.keys(buildAggregateRows(runInfo)[0]), buildAggregateRows(runInfo));
  writeSheet(workbook, "g12_g13_nonhdr_data", Object.keys(buildLongRows(tsvRows, figureReadyPath, runInfo)[0]), buildLongRows(tsvRows, figureReadyPath, runInfo));
  writeSheet(workbook, "per_sample_summary", Object.keys(tsvRows[0]), tsvRows);
  writeSheet(workbook, "removed_hdr_summary", Object.keys(buildRemovedHdrRows(runInfo)[0]), buildRemovedHdrRows(runInfo), { headerFill: "#7F6000" });
  writeSheet(workbook, "qc_checks", Object.keys(buildQcRows(runInfo)[0]), buildQcRows(runInfo), { headerFill: "#38761D" });
  writeSheet(workbook, "provenance", Object.keys(buildProvenanceRows(runInfo)[0]), buildProvenanceRows(runInfo), { headerFill: "#134F5C" });

  setUsefulWidths(workbook);
  await compactVerify(workbook);

  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(workbookPath);

  runInfo.outputs.workbook_xlsx = workbookPath;
  runInfo.checksums_sha256.outputs.workbook_xlsx = await sha256(workbookPath);
  await fs.writeFile(runInfoPath, JSON.stringify(runInfo, null, 2) + "\n", "utf8");
  const manifestPath = await writeManifest(packageDir, runInfo, workbookPath);

  console.log(JSON.stringify({
    workbook_path: workbookPath,
    manifest_path: manifestPath,
    workbook_sha256: runInfo.checksums_sha256.outputs.workbook_xlsx,
  }, null, 2));
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
