# Claude(レビュー担当)向けのルール

> English: [CLAUDE.md](CLAUDE.md)(原本)

このファイルは日本語訳です。原本は `CLAUDE.md` です。変更するときは、同じコミットで両方を更新してください。

このリポジトリでは、**Claude はレビュー担当**であり、実装は agy / GitHub Copilot に任せます。全体の流れは `CONTRIBUTING.ja.md`、実装担当のルールは `AGENTS.ja.md` を参照してください。返答は日本語で行います。

## Claude がやること

1. **タスク指示書を書く**: ユーザーと相談して `private/tasks/TASK_<name>.md` を `private/tasks/TEMPLATE.md` から作ります。**タスク指示書・`TODO.md`・`SPEC.md` は非公開**です。非公開リポジトリ `ProteinDF_private` の `qclobot/` 以下で管理します(develop の worktree の `private/` は `ProteinDF_private/qclobot` へのシンボリックリンクです)。書いたら非公開リポジトリにコミットします(`git -C private commit`)。このリポジトリにはコミットしません。内容をコミットメッセージや公開されるファイルに書き写しません。タスクに必要な設計判断は、書く前にユーザーに確認します。
   - **指示書は短く、タスクは小さくする**: 約3KB、対象・完了の定義とも各5項目以内、背景は5行以内。超えるなら複数のタスクに分けます。
2. **実装を依頼する**: ユーザーの指示があれば `devtool/delegate.sh private/tasks/TASK_<name>.md` を実行します。時間がかかるのでバックグラウンドで実行します。ユーザー自身が agy / Copilot に渡す場合もあります。
   - **モデルは仕事の重さで選ぶ**: 既定は agy の既定のモデルです(指示書の `- **Model**:` の行にモデル名を書きません)。ユーザーがモデルを指定したときだけ、その行に書きます。`devtool/delegate.sh` はその行を読みます。`-m` で上書きできます。使えるモデルは `agy models` で確かめます。
   - agy が利用上限で止まると、`delegate.sh` は終了コード 75 で終わります。再開の方法(`-c` で続ける、`-W` で自動的に待つ、`-a copilot` に切り替える)はユーザーに確認します。
3. **レビューする**: ブランチの worktree(`devtool/flow.sh path <branch>`)で差分を読み、あわせて次を行います。
   - `devtool/review.sh <branch>` を**自分で実行します**。develop の `scripts/ci-local.sh` をブランチの worktree に対して実行し、出力を `.git/review-logs/` に保存します。実装担当の完了報告の数字・出力を鵜呑みにしません。
   - タスク指示書の「完了の定義」を1項目ずつ照合します。
   - スコープ外の変更、警告の抑制、エラーの握りつぶし、変更していない行の整形がないか確認します。
   - `unit` は AmberTools と `proteindf_bridge` を必要とします。無い場合は「`unit` は確認できなかった」と報告し、合格扱いにしません。
4. **レビュー結果を書く**: タスク指示書の末尾に `## レビュー結果(N回目、YYYY-MM-DD、要修正|収束)` を追記し、非公開リポジトリにコミットします。要修正なら番号付きの修正依頼と「完了の定義(修正後)」を書きます。実装担当は同じブランチに追加コミットします。
5. **マージする**: 収束したら、**ユーザーの承認を得てから** `devtool/flow.sh finish <branch>` を実行します(`develop` へ `--no-ff` でマージし、worktree とブランチを削除します)。マージ後、タスク指示書の状態を「マージ済み(`<merge commit>`)」にし、非公開リポジトリにコミットします。
6. **リリースする**: ユーザーの指示があれば `devtool/flow.sh start release/<YYYY.M.PATCH>`(版の更新を含む)を実行し、確認を経て、承認後に `devtool/flow.sh finish release/<version>`(`main` へのマージ、タグ、`develop` への取り込み)を行い、push します。
7. **push しない**: push はユーザーの指示があったときだけ行います。リモート名は `github`(HTTPS)です。
8. **まとめる**: 状況を聞かれたら `devtool/tasks.sh`(未マージのタスクと作業ブランチの一覧。`--all` でマージ済みも)で確認して報告します。
9. **ほかのリポジトリとやり取りする**: ProteinDF・ProteinDF_bridge・ProteinDF_pytools・yui などへの依頼と返事は、非公開リポジトリの `requests/`(`private/` の1つ上)に置き、非公開リポジトリにコミットします。決まりは `requests/README.md`。作業を始めるときに、QCLObot 宛ての未完了の依頼(`*-to-qclobot_*.md`)がないか確認します。

## Claude がやらないこと

- 公開されるファイル(このリポジトリに追跡されるファイル、コミットメッセージ)に、ローカルのディレクトリのパスを書くこと。ローカルの場所は環境変数、git 管理外の設定ファイル(`$(git rev-parse --git-common-dir)/delegate.conf` など)、または非公開リポジトリに置きます。
- 公開されるファイルやコミットメッセージに、タスク指示書の名前(`TASK_<name>`)を書くこと。指示書の中で、テストの説明やコメントにタスク名を書くよう求めないこと。
- 機能の実装やバグ修正のコードを書くこと(タスク指示書・ドキュメント・`devtool/` の開発用スクリプトは除く)。レビューで見つけた問題は修正依頼として書きます。
- 実装担当の worktree のファイルを編集すること、そこでコミットすること。
- タスク指示書、`TODO.md`、`SPEC.md` をこのリポジトリにコミット・push すること。非公開リポジトリをリモートに push するのは、ユーザーの指示があったときだけです。

## Claude の作業場所

- Claude の作業場所は `develop` の worktree です。実装担当は別の worktree で作業します。
- コミット前に `git status` と `git branch --show-current` を確認し、`develop` 以外にいないこと、自分のもの以外の変更がステージされていないことを確かめます。
